%% Sistema de 15 Barras - Tabu Search

% std = 0 para 100 vezes!
% Algoritmo: Heurística obtém solução inicial, quando
% há estagnação entra a Busca Tabu.
% de início há a necessidade da solução ser factível,
% depois essa restrição é removida.

%% Dados dos relés

% Plug Settings
% PS_otim = [1.50 1.00 2.00 1.00 2.00 2.00 2.00 1.50 2.00 1.50 1.50 1.50 2.00 1.00 1.00 1.50 2.00 1.00 2.00 1.50 0.50 1.50 1.00 1.50 2.00 1.50 2.00 2.50 1.50 2.00 2.00 1.50 2.50 2.50 2.00 2.00 2.50 2.50 2.50 2.50 2.50 1.50]; % otimização Kida
PS_disp = zeros(count_MAX,1);

for i = 1:count_MAX
    PS_disp(i,1) = 1.5 + 0.5*(i-1);
end

% Geração de Ponto Inicial Aleatório
PS = zeros(N_Reles,1);

rng('shuffle');
for i = 1:N_Reles
    %     PS(i,1) = PS_disp(randi(5),1);
    %     PS_inicial(i,vezes) = PS(i,1);
    PS(i,1) = PS_inicial(i,vezes);        % Para PS específicos (teste)
end


% Flags
PL_count = 0;

CTI = 0.2;
RTC = [300/5, 200/5, 200/5, 300/5, 200/5, 400/5];

% Ordem dos relês (prim-p / ret-q)
p = [1; 2; 3; 4; 5; 6];
q = [5; 4; 1; 6; 3; 2];


% Faltas modo normal
FC_pr = [1978.9 1525.7 1683.9 1815.4 1499.66 1766.3];
FC_bc = [175 545 617.22 466.17 384 145.34];

% Faltas modo transiente
FC_pr_Tr = [2075 1621.7 1779.6 1911.5 1588.5 1855.4];
FC_bc_Tr = [400.7 700.64 760.17 622.65 558.13 380.7];

Icc_K1 = [FC_pr; FC_bc];    % modo normal
Icc_K2 = [FC_pr_Tr; FC_bc_Tr];    % modo transiente

B1 = size(Icc_K1);
B2 = size(Icc_K2);
f = zeros(1,N_Reles);

% Curva IEC Inversa
Ai = 0.14;
N = 0.02;

options = optimoptions(@linprog,'Display','none');

%% Vetores Tabu

% Critérios de Frequência
moves_usage = zeros(count_MAX, N_Reles);      % Valores de PS mais utilizados (Diversificação)
%more_frequent_moves = zeros(1,N_Reles);
tabu_short_mem = zeros(count_MAX, N_Reles);           % valores de PS proibidos, dimensão (PS_disp x Relês)

% Melhor movimento
exist_flag = 0;                         % Existe melhor movimento disponível?
best_move = ones(N_Reles,1) * 1.5;
best_B1 = size(B1);
best_B2 = size(B2);
best_fo = 10000;                        % melhor resultado
best_fo_it = 10000;                     % para definir parada
pos_tabu_best = 0;                      % posição do melhor movimento

% (Inativos)
pos_tabu_worst = 0;                     % posição do pior movimento
worst_move = ones(N_Reles,1) * 1.5;
worst_fo = 0.1;

% Bad Move (Diversificação apenas!)
exist_best_move = zeros(N_Reles,1);          % Indica se há um melhor movimento disponível
bad_move = ones(N_Reles,1) * 1.5;            % Movimento que menos piora a F.O
bad_fo = 1000;
bad_fo_it = 0;
bad_B1 = size(B1);
bad_B2 = size(B2);

% Iterações
%it_max = 20;
it_results = zeros(200,2);              % F.O
it_lim_div = 0;                         % Iteração máx. estágio de Div.
it_lim_it = 0;                          % Iteração máx. estágio de Int.
Tabu_Search_it = 0;                             % Iteração na qual foi pra Diversificação
Heuristic_it = 0;


% Duração do banimento
tabu_tenure_best = 0;                   % banimento do melhor movimento
tabu_tenure_worst = 0;                  % banimento do pior movimento (inativo)

tabu_state = 0;                         % Estágio do algoritmo:
%   0: Intensificação (busca)
%   1: Diversificação

%% Resultados

%% Loop

for it = 1:N_it_MAX
    
    switch tabu_state
        
        case 0
            %% Intensificação
            for R = 1:N_Reles
                for count = 1:count_MAX
                    if(tabu_short_mem(count,R) == 0)    % Somente se o vetor não é tabu
                        % que deve-se calcular o PL!
                        % Atualização do PS (incremento)
                        if((R ~= 1) || (count ~= 1))
                            PS(R,1) = PS_disp(count,1);
                        end
                        
                        % Matrizes Kij
                        for j = 1:N_Reles
                            % Matriz faltas normais
                            B1(1,j) = Ai/((FC_pr(j)/(PS(j)*RTC(j)))^N-1);
                            B1(2,q(j)) = Ai/((FC_bc(j)/(PS(q(j))*RTC(q(j))))^N-1);
                            % Matriz transiente
                            B2(1,j) = Ai/((FC_pr_Tr(j)/(PS(j)*RTC(j)))^N-1);
                            B2(2,q(j)) = Ai/((FC_bc_Tr(j)/(PS(q(j))*RTC(q(j))))^N-1);
                        end
                        
                        Kij_geral = zeros(N_Restricoes,N_Reles);
                        
                        for k = 1:6
                            Kij_geral(k,p(k)) = B1(1,p(k));
                            Kij_geral(k,q(k)) = -B1(2,q(k));
                            Kij_geral(k+6,p(k)) = B2(1,p(k));
                            Kij_geral(k+6,q(k)) = -B2(2,q(k));
                        end
                        
                        %% Sistema no formato do Linprog
                        
                        % F.O - Tempo mínimo de atuação dos Relés (somente Primários)
                        
                        for i = 1:N_Reles
                            % Somente um tipo de falta é levado em conta no
                            % cálculo da FO!
                            f(1,i) = B1(1,i);
                        end
                        
                        % Valores Mín/Máx de Dial dos Relés R1-R6
                        lb = ones(1,N_Reles);
                        lb = lb * 0.1;
                        ub = ones(1,N_Reles);
                        ub = ub * 1.1;
                        
                        % Lado direito da Desigualdade
                        b = ones(N_Restricoes,1);
                        b = b * -CTI;
                        
                        A = Kij_geral;
                        
                        % NÃO HÁ EXCEÇÕES NO SISTEMA DE 3 BARRAS!
                        
                        %% Resolução via PL
                        [x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);
                        PL_count = PL_count +1;
                        
                        %% Melhor movimento e atualização das Listas Tabu
                        if((exitflag == 1) || (it > 2))
                            if(fval <= best_fo)
                                best_B1 = B1;
                                best_B2 = B2;
                                best_f = f;
                                best_x = x;
                                best_move(R,1) = PS(R,1);
                                best_fo = fval;
                                pos_tabu_best = count;
                                exist_best_move(R,1) = -1;
                                moves_usage(count,R) = moves_usage(count,R) +1;
                            else
                                if(fval <= bad_fo)
                                    % Caso piore a F.O mas seja a melhor opção
                                    % disponível, armazenar o valor
                                    bad_B1 = B1;
                                    bad_B2 = B2;
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
                            % Voltar para última configuração factível.
                            PS = best_move;
                        end
                        
                        if(count == count_MAX)
                            % Atualização da Lista Tabu (decremeneto)
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
                            
                            % Caso esteja estagnado, muda estratégia de
                            % busca
                            if(Tabu_Search_it ~= 0)
                                if((it <= Tabu_Search_it + 3) && (Heuristic_it == 0))
                                    tabu_tenure_best = 1;
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
                            % Atualização da Lista Tabu (decremeneto)
                            for i = 1:count_MAX
                                if(tabu_short_mem(i,R) > 0)
                                    tabu_short_mem(i,R) = tabu_short_mem(i,R) - 1;
                                end
                            end
                        end
                    end
                end
            end
            
            %% Diversificação
        case 1
            % Selecionar movimento mais frequente.
            %                         [M,Ind] = max(moves_usage,[],1);
            %                         for R = 1:N_Reles
            %                             PS(R,1) = PS_disp(Ind(1,R),1);
            %                         end
            %                         clear M Ind;
            tabu_state = 0;
            
    end
    % Progresso das iterações
    if(tabu_state == 0)
        it_results(it,1) = best_fo;
        it_results(it,2) = bad_fo_it;
        fprintf('Run = %d, it = %d, FO = %f, Tabu_Search_it = %d, Heuristic_it = %d \n',vezes, it, best_fo, Tabu_Search_it, Heuristic_it);
    end
    
    
    
    if(it >= 3 + Tabu_Search_it)
        penult_result = it_results((it-2),1);
        last_result = it_results((it-1),1);
        current_result = it_results((it),1);
        
        % Se a Busca está estagnada por 3 iterações, inicia Busca Tabu
        if(current_result == last_result)
            if(current_result == penult_result)
                % Caso já tenha executado Busca Tabu 1x e estagnou
                % novamente, volta para Heurística
                if(Tabu_Search_it == 0)
                    Tabu_Search_it = it;        % Ativa a Busca Tabu
                else if(Heuristic_it == 0)    % Desativa Busca Tabu, Volta pra Heurística
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


%% Recalcular Função Objetivo com melhores valores

% Melhores valores
PS = best_move;
B1 = best_B1;
B2 = best_B2;
f = best_f;

for j = 1:6
    % Matriz faltas normais
    B1(1,j) = Ai/((FC_pr(j)/(PS(j)*RTC(j)))^N-1);
    B1(2,q(j)) = Ai/((FC_bc(j)/(PS(q(j))*RTC(q(j))))^N-1);
    % Matriz transiente
    B2(1,j) = Ai/((FC_pr_Tr(j)/(PS(j)*RTC(j)))^N-1);
    B2(2,q(j)) = Ai/((FC_bc_Tr(j)/(PS(q(j))*RTC(q(j))))^N-1);
end

for k = 1:6
    Kij_geral(k,p(k)) = B1(1,p(k));
    Kij_geral(k,q(k)) = -B1(2,q(k));
    Kij_geral(k+6,p(k)) = B2(1,p(k));
    Kij_geral(k+6,q(k)) = -B2(2,q(k));
end
%% Sistema no formato do Linprog

% F.O - Tempo mínimo de atuação dos Relés

A = Kij_geral;

% NÃO HÁ EXCEÇÕES NO SISTEMA DE 3 BARRAS!

[x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);



%% Tempos de Atuação dos Relês

T_pri1 = zeros(1,N_Reles);    % Pré-alocação
T_ret1 = zeros(1,N_Reles);
T_pri2 = zeros(1,N_Reles);
T_ret2 = zeros(1,N_Reles);

for i = 1:N_Reles
    T_pri1(1,i) = x(i,1) * B1(1,i);
    T_ret1(1,i) = x(i,1) * B1(2,i);
    T_pri2(1,i) = x(i,1) * B2(1,i);
    T_ret2(1,i) = x(i,1) * B2(2,i);
end

T_coord = ones(N_Restricoes,8);     % 2 Configurações: Normal e Transiente.


violation = zeros(N_Restricoes,3);
index = 1;

for i = 1:N_Restricoes/2
    T_coord(i,1) = p(i);
    T_coord(i,2) = q(i);
    T_coord(i,3) = T_pri1(p(i));
    T_coord(i,4) = T_ret1(q(i));
    T_coord(i,5) = T_pri2(p(i));
    T_coord(i,6) = T_ret2(q(i));
    T_coord(i,7) = T_coord(i,4) - T_coord(i,3);
    T_coord(i,8) = T_coord(i,6) - T_coord(i,5);
end

t = toc;