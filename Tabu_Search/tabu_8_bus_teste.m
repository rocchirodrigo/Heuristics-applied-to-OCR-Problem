%% Sistema de 8 Barras - Tabu Search

% std = XXX para 100 vezes!
% Algoritmo: Heur�stica obt�m solu��o inicial, quando
% h� estagna��o entra a Busca Tabu.
% de in�cio h� a necessidade da solu��o ser fact�vel,
% depois essa restri��o � removida.

%% Dados dos rel�s

% Plug Settings
% PS_otim = [1.50 1.00 2.00 1.00 2.00 2.00 2.00 1.50 2.00 1.50 1.50 1.50 2.00 1.00 1.00 1.50 2.00 1.00 2.00 1.50 0.50 1.50 1.00 1.50 2.00 1.50 2.00 2.50 1.50 2.00 2.00 1.50 2.50 2.50 2.00 2.00 2.50 2.50 2.50 2.50 2.50 1.50]; % otimiza��o Kida

PS_disp = [0.5 0.6 0.8 1.0 1.5 2.0 2.5]';

% Gera��o de Ponto Inicial Aleat�rio
PS = zeros(N_Reles,1);

rng('shuffle');
for i = 1:N_Reles
    PS(i,1) = PS_disp(randi(count_MAX),1);
    PS_inicial(i,vezes) = PS(i,1);
%     PS(i,1) = PS_inicial(i,vezes);        % Para PS espec�ficos (teste)
end

% Flags
PL_count = 0;

CTI = 0.3;
RTC = [240, 240, 160, 240, 240, 240, 160, 240, 160, 240, 240, 240, 240, 160];

% Ordem dos rel�s (prim-p / ret-q)
p = [1;2;2;3;4;5;6;6;7;7;8;8;9;10;11;12;12;13;14;14];

q = [6;1;7;2;3;4;5;14;5;13;7;9;10;11;12;13;14;8;1;9];

% N�veis de Corrente de Falta.
% Sem exce��es!

FC_pr = [3232, 5924, 3556, 3783, 2401, 6109, 5223, 6093, 2484, 3883, 3707, 5899, 2991, 5199];

FC_bc = [996, 3556, 2244, 2401, 1197, 3232, 1890, 2991, 1165, 2484, 2344, 3707, 987, 1874];

Icc_K1 = [FC_pr; FC_bc];

B1 = size(Icc_K1);
f = zeros(1,N_Reles);

% Curva IEC Inversa
Ai = 0.14;
N = 0.02;

options = optimoptions(@linprog,'Display','none');

%% Vetores Tabu

% Crit�rios de Frequ�ncia
moves_usage = zeros(count_MAX, N_Reles);      % Valores de PS mais utilizados (Diversifica��o)
%more_frequent_moves = zeros(1,N_Reles);
tabu_short_mem = zeros(count_MAX, N_Reles);           % valores de PS proibidos, dimens�o (PS_disp x Rel�s)

% Melhor movimento
exist_flag = 0;                         % Existe melhor movimento dispon�vel?
best_move = ones(N_Reles,1) * 0.5;
best_B1 = size(B1);
best_fo = 10000;                        % melhor resultado
best_fo_it = 10000;                     % para definir parada
pos_tabu_best = 0;                      % posi��o do melhor movimento

% (Inativos)
pos_tabu_worst = 0;                     % posi��o do pior movimento
worst_move = ones(N_Reles,1) * 0.5;
worst_fo = 0.1;

% Bad Move (Diversifica��o apenas!)
exist_best_move = zeros(N_Reles,1);          % Indica se h� um melhor movimento dispon�vel
bad_move = ones(N_Reles,1) * 0.5;            % Movimento que menos piora a F.O
bad_fo = 1000;
bad_fo_it = 0;
bad_B1 = size(B1);

% Itera��es
%it_max = 20;
it_results = zeros(200,2);              % F.O
it_lim_div = 0;                         % Itera��o m�x. est�gio de Div.
it_lim_it = 0;                          % Itera��o m�x. est�gio de Int.
Tabu_Search_it = 0;                             % Itera��o na qual foi pra Diversifica��o
Heuristic_it = 0;


% Dura��o do banimento
tabu_tenure_best = 0;                   % banimento do melhor movimento
tabu_tenure_worst = 0;                  % banimento do pior movimento (inativo)

tabu_state = 0;                         % Est�gio do algoritmo:
%   0: Intensifica��o (busca)
%   1: Diversifica��o

%% Resultados

%% Loop

for it = 1:N_it_MAX
    
    switch tabu_state
        
        case 0
            %% Intensifica��o
            for R = 1:N_Reles
                for count = 1:count_MAX
                    if(tabu_short_mem(count,R) == 0)    % Somente se o vetor n�o � tabu
                        % que deve-se calcular o PL!
                        % Atualiza��o do PS (incremento)
                        if((R ~= 1) || (count ~= 1))
                            PS(R,1) = PS_disp(count,1);
                        end
                        
                        % Matrizes Kij
                        for j = 1:N_Reles
                            B1(1,j) = Ai/((FC_pr(j)/(PS(j)*RTC(j)))^N-1);  % Rel� prim�rio
                            B1(2,j) = Ai/((FC_bc(j)/(PS(j)*RTC(j)))^N-1);  % Rel� retaguarda
                        end
                        
                        Kij_geral = zeros(N_Restricoes,N_Reles);
                        
                        for k = 1:N_Restricoes
                            Kij_geral(k,p(k)) = B1(1,p(k));
                            Kij_geral(k,q(k)) = -B1(2,q(k));
                        end
                        
                        %% Sistema no formato do Linprog
                        
                        % F.O - Tempo m�nimo de atua��o dos Rel�s (somente Prim�rios)
                        
                        for i = 1:N_Reles
                            f(1,i) = B1(1,i);
                        end
                        
                        % Valores M�n/M�x de Dial dos Rel�s R1-R6
                        lb = ones(1,N_Reles);
                        lb = lb * 0.1;
                        ub = ones(1,N_Reles);
                        ub = ub * 1.1;
                        
                        % Lado direito da Desigualdade
                        b = ones(N_Restricoes,1);
                        b = b * -CTI;
                        
                        A = Kij_geral;
                        
                        % Exce��es: N�O H� EXCE��ES!
                        
                        %% Resolu��o via PL
                        [x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);
                        PL_count = PL_count +1;
                        
                        %% Melhor movimento e atualiza��o das Listas Tabu
                        if((exitflag == 1) || (it > 2))
                            if(fval <= best_fo)
                                best_B1 = B1;
                                best_f = f;
                                best_x = x;
                                best_move(R,1) = PS(R,1);
                                best_fo = fval;
                                pos_tabu_best = count;
                                exist_best_move(R,1) = -1;
                                moves_usage(count,R) = moves_usage(count,R) +1;
                            else
                                if(fval <= bad_fo)
                                    % Caso piore a F.O mas seja a melhor op��o
                                    % dispon�vel, armazenar o valor
                                    bad_B1 = B1;
                                    bad_f = f;
                                    bad_x = x;
                                    bad_move(R,1) = PS(R,1);
                                    bad_fo = fval;
                                end
                                if(fval > worst_fo)
                                    worst_fo = fval;
                                    worst_move(R,1) = PS(R,1);
                                    pos_tabu_worst = count;
                                end
                            end
                        else
                            % Voltar para �ltima configura��o fact�vel.
                            PS = best_move;
                        end
                        
                        if(count == count_MAX)
                            % Atualiza��o da Lista Tabu (decremeneto)
                            for i = 1:count_MAX
                                if(tabu_short_mem(i,R) > 0)
                                    tabu_short_mem(i,R) = tabu_short_mem(i,R) - 1;
                                end
                            end
                            
                            PS(R,1) = best_move(R,1);
                            
                            % Melhor movimento
                            if(pos_tabu_best ~=0)
                                tabu_short_mem(pos_tabu_best,R) = tabu_tenure_best;
                            end
                            
                            % Pior movimento
                            if(pos_tabu_worst ~=0)
                                tabu_short_mem(pos_tabu_worst,R) = tabu_tenure_worst;
                            end
                            
                            % Caso esteja estagnado, muda estrat�gia de
                            % busca
                            if(Tabu_Search_it ~= 0)
                                if((it <= Tabu_Search_it + 3) && (Heuristic_it == 0))
                                    tabu_tenure_best = 5;
                                    if(exist_best_move(R,1) == 0)
                                        PS(R,1) = bad_move(R,1);
                                        bad_fo_it = bad_fo;
                                        bad_fo = 1000;
                                        exist_flag = exist_flag + 1;
                                    end
                                else tabu_tenure_best = 0;
                                end
                            end
                            
                            % Zera o vetor
                            for i = 1:N_Reles
                                exist_best_move(i,1) = 0;
                            end
                        end
                    else
                        if(count == count_MAX)
                            % Atualiza��o da Lista Tabu (decremeneto)
                            for i = 1:count_MAX
                                if(tabu_short_mem(i,R) > 0)
                                    tabu_short_mem(i,R) = tabu_short_mem(i,R) - 1;
                                end
                            end
                        end
                    end
                end
            end
            
            %% Diversifica��o
        case 1
            % Selecionar movimento mais frequente.
%                         [M,Ind] = max(moves_usage,[],1);
%                         for R = 1:N_Reles
%                             PS(R,1) = PS_disp(Ind(1,R),1);
%                         end
%                         clear M Ind;
                        tabu_state = 0;
            
    end
    % Progresso das itera��es
    if(tabu_state == 0)
        it_results(it,1) = best_fo;
        it_results(it,2) = bad_fo_it;
        fprintf('Run = %d, it = %d, FO = %f, Tabu_Search_it = %d, Heuristic_it = %d \n',vezes, it, best_fo, Tabu_Search_it, Heuristic_it);
    end
    
    
    % Se a Busca est� estagnada por 3 itera��es, inicia Busca Tabu
    if(it >= 3 + Tabu_Search_it)
        penult_result = it_results((it-2),1);
        last_result = it_results((it-1),1);
        current_result = it_results((it),1);
        if(current_result == last_result)
            if(current_result == penult_result)
                % Caso j� tenha executado Busca Tabu 1x e estagnou
                % novamente, volta para Heur�stica
                if(Tabu_Search_it == 0)
                    Tabu_Search_it = it;        % Ativa a Busca Tabu
                else if(Heuristic_it == 0)    % Desativa Busca Tabu, Volta pra Heur�stica
                        Tabu_Search_it = 0;
                        Heuristic_it = it;
                    else if(it >= 3 + Heuristic_it)
                            tabu_state = 1;
                        end
                    end
                end
            end
        end
    end
end


%% Recalcular Fun��o Objetivo com melhores valores

% Melhores valores
PS = best_move;
B1 = best_B1;
f = best_f;

for k = 1:N_Restricoes
    Kij_geral(k,p(k)) = B1(1,p(k));
    Kij_geral(k,q(k)) = -B1(2,q(k));
end

%% Sistema no formato do Linprog

% F.O - Tempo m�nimo de atua��o dos Rel�s

A = Kij_geral;

% Exce��es: N�O H� EXCE��ES!

[x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);

%% Tempos de Atua��o dos Rel�s

T_pri = zeros(N_Restricoes,1);    % Pr�-aloca��o
T_ret = zeros(N_Restricoes,1);

for i = 1:N_Reles
    T_pri(i,1) = x(i,1) * B1(1,i);
    T_ret(i,1) = x(i,1) * B1(2,i);
end

T_coord = ones(N_Restricoes,5);     % Configura��o �nica!
violation = zeros(N_Restricoes,3);
index = 1;

for i = 1:N_Restricoes
    T_coord(i,1) = p(i);
    T_coord(i,2) = q(i);
    T_coord(i,3) = T_pri(p(i));
    T_coord(i,4) = T_ret(q(i));
    T_coord(i,5) = T_coord(i,4) - T_coord(i,3);
    
    if(T_coord(i,5) < CTI)  % Viola��es na Coordena��o
        violation(index,1) = p(i);
        violation(index,2) = q(i);
        violation(index,3) = T_coord(i,5);
        index = index+1;
    end
end

t = toc;