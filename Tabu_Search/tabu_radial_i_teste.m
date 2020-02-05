%% Sistema I - Tabu Search

%% COLOCAR DADOS DO SISTEMA I!!

% Algoritmo: Heurística obtém solução inicial, quando
% há estagnação entra a Busca Tabu.
% de início há a necessidade da solução ser factível,
% depois essa restrição é removida.

%% Dados dos relés

% Plug Settings
PS_disp = zeros(count_MAX,1);

for i = 1:count_MAX
    PS_disp(i,1) = 2.5 + 0.75*(i-1);
end

% Geração de Ponto Inicial Aleatório
PS = zeros(N_Reles,1);

rng('shuffle');
for i = 1:N_Reles
    PS(i,1) = PS_disp(randi(count_MAX),1);
    PS_inicial(i,vezes) = PS(i,1);
%         PS(i,1) = PS_inicial(i,vezes);        % Para PS específicos (teste)
end

% Flags
PL_count = 0;

CTI = 0.4;
RTC = [300/5, 300/5, 100/5, 200/5, 100/5];

% Ordem dos relês (prim-p / ret-q)
p = [1; 2; 3; 4; 5];
q = [1; 1; 1; 2; 3];

% Odem dos relês primários na matriz
r = [1; 2; 3; 4];
s = [2; 3; 4; 5];

% Odem dos relês retaguarda na matriz
u = [1; 1; 2; 3];

% Faltas máx-mín
FC_max = [3115.0 2010.7 2010.7 1512.5 878.4];
FC_min = [1741.3 1309.9 760.7 500.3 325.1];

Icc_K1 = FC_max;    % Icc máx
Icc_K2 = FC_min;    % Icc mín

% Corrente de Carga
I_load = [199.5 130.8 68.7 100.7 50.0];

% Limites PS - Carga
FC = 1.50;
lim_inf = ones(N_Reles,1);
lim_sup = ones(N_Reles,1);

% Restrição do valor de PS
for i = 1:N_Reles
    lim_inf(i,1) = (FC*I_load(i))/RTC(i);
    lim_sup(i,1) = FC_min(i)/RTC(i);
end

% Valores mínimos permitidos pelo Relé
PS_min = interp1(PS_disp, PS_disp, lim_inf,'next');
PS_max = interp1(PS_disp, PS_disp, lim_sup,'nearest','extrap');
rest_carga = ones(count_MAX,N_Reles);
rest_carga2 = ones(count_MAX,N_Reles);

% Ajusta o Ponto Inicial para valores permitidos pelo Relé
for i = 1:N_Reles
    if(PS(i,1) < PS_min(i,1))
        PS(i,1) = PS_min(i,1);
    else if(PS(i,1) > PS_max(i,1))
            PS(i,1) = PS_max(i,1);
        end
    end
end

% Restrição para evitar cálculo de PLs quando PS está fora do intervalo
% permitido
for i = 1:count_MAX
    for j = 1:N_Reles
        rest_carga(i,j) = PS_disp(i,1);
        rest_carga2(i,j) = PS_disp(i,1);
        if((rest_carga(i,j)) < lim_inf(j,1) || (rest_carga(i,j)) > lim_sup(j,1))
            rest_carga(i,j) = 1;    % Fora da restrição
        else rest_carga(i,j) = 0;   % Cumpre a restrição
        end
    end
end

B1 = size(Icc_K1);
B2 = size(Icc_K2);
f = zeros(1,N_Reles);

% Curva IEC Inversa
Ai = [13.50 0.14 0.14 13.50 80.00];
N = [1.00 0.02 0.02 1.00 2.00];

options = optimoptions(@linprog,'Display','none');

%% Vetores Tabu

% Critérios de Frequência
moves_usage = zeros(count_MAX, N_Reles);      % Valores de PS mais utilizados (Diversificação)
%more_frequent_moves = zeros(1,42);
tabu_short_mem = zeros(count_MAX, N_Reles);           % valores de PS proibidos, dimensão (PS_disp x Relês)

% Melhor movimento
exist_flag = 0;                         % Existe melhor movimento disponível?

% Melhor movimento inicializado no menor valor de PS permitido.
best_move = zeros(N_Reles,1);
for i = 1:N_Reles
    best_move(i,1) = PS_min(i,1);
end

best_fo = 10000;                        % melhor resultado
best_fo_it = 10000;                     % para definir parada
pos_tabu_best = 0;                      % posição do melhor movimento

% Bad Move (Diversificação apenas!)
exist_best_move = zeros(N_Reles,1);          % Indica se há um melhor movimento disponível
bad_move = PS_min;                           % Movimento que menos piora a F.O
bad_fo = 1000;
bad_fo_it = 0;
bad_B1 = size(B1);
bad_B2 = size(B2);

% Iterações
it_results = zeros(200,2);              % F.O
it_lim_div = 0;                         % Iteração máx. estágio de Div.
it_lim_it = 0;                          % Iteração máx. estágio de Int.
Tabu_Search_it = 0;                             % Iteração na qual foi pra Diversificação
Heuristic_it = 0;


% Duração do banimento
tabu_tenure_best = 0;                   % banimento do melhor movimento
tabu_tenure_worst = 0;                  % banimento do pior movimento (inativo)

tabu_state = 0;                         % Estágio do algoritmo:
%   0: Heurística
%   1: Busca Tabu

%% Resultados

%% Loop

for it = 1:N_it_MAX
    
    switch tabu_state
        case 0
            %% Intensificação
            for R = 1:N_Reles
                for count = 1:count_MAX
                    if((tabu_short_mem(count,R) == 0) && rest_carga(count,R) == 0)
                        % Somente se o vetor não é tabu e não está restrito
                        % por valor limite do PS que deve-se calcular o PL!
                        
                        % Atualização do PS (incremento)
                        if((R ~= 1) || (count ~= 1))
                            PS(R,1) = PS_disp(count,1);
                        end
                        
                        % Matrizes Kij
                        for j = 1:N_Reles
                            % Matriz faltas max
                            B1(1,j) = Ai(j)/((FC_max(j)/(PS(j)*RTC(j)))^N(j)-1);
                            B1(2,j) = Ai(q(j))/((FC_max(j)/(PS(q(j))*RTC(q(j))))^N(q(j))-1);
                            % Matriz faltas min
                            B2(1,j) = Ai(j)/((FC_min(j)/(PS(j)*RTC(j)))^N(j)-1);
                            B2(2,j) = Ai(q(j))/((FC_min(j)/(PS(q(j))*RTC(q(j))))^N(q(j))-1);
                        end
                        
                        Kij_geral = zeros(N_Restricoes,N_Reles);
                        
                        for k = 1:4
                            Kij_geral(k,s(k)) = B1(1,s(k)); % prim. max
                            Kij_geral(k,u(k)) = -B1(2,s(k)); % ret. max
                            Kij_geral(4+k,s(k)) = B2(1,s(k)); % prim. min
                            Kij_geral(4+k,u(k)) = -B2(2,s(k)); % ret. min
                        end
                        
                        for k = 1:5
                            % Tmin
                            Kij_geral(8+k,k) = -B1(1,k);
                            Kij_geral(13+k,k) = -B2(1,k);
                            % Tmax
                            Kij_geral(18+k,k) = B1(1,k);
                            Kij_geral(23+k,k) = B2(1,k);
                        end
                        
                        %% Sistema no formato do Linprog
                        
                        % F.O - Tempo mínimo de atuação dos Relés (conta
                        % valores de curto mínimo e máximos!)
                        
                        for i = 1:N_Reles
                            f(1,i) = B1(1,i) + B2(1,i);
                        end
                        
                        % Valores Mín/Máx de Dial dos Relés R1-R5
                        lb = ones(1,N_Reles) * 0.1;
                        ub = ones(1,N_Reles) * 10;
                        
                        % Lado direito da Desigualdade
                        b = ones(N_Restricoes,1);
                        
                        % CTI
                        for i = 1:8
                            b(i,1) = b(i,1) * -0.4;
                        end
                        
                        %T min
                        for i = 1:10
                            b(i+8,1) = b(i+8,1) * -0.05;
                        end
                        
                        % Tmax
                        for i = 1:10
                            b(i+18,1) = b(i+18,1) * 2;
                        end
                        A = Kij_geral;
                        
                        %% Resolução via PL
                        [x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);
                        PL_count = PL_count +1;
                        
                        %% Melhor movimento e atualização das Listas Tabu
                        if(exitflag == 1)
                            if(fval <= best_fo)
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
                                    bad_f = f;
                                    bad_x = x;
                                    bad_move(R,1) = PS(R,1);
                                    bad_fo = fval;
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
                            
%                             PS(R,1) = best_move(R,1);
                            
                            % Melhor movimento
                            if(pos_tabu_best ~=0)
                                tabu_short_mem(pos_tabu_best,R) = tabu_tenure_best;
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

for j = 1:N_Reles
    % Matriz faltas max
    B1(1,j) = Ai(j)/((FC_max(j)/(PS(j)*RTC(j)))^N(j)-1);
    B1(2,j) = Ai(q(j))/((FC_max(j)/(PS(q(j))*RTC(q(j))))^N(q(j))-1);
    % Matriz faltas min
    B2(1,j) = Ai(j)/((FC_min(j)/(PS(j)*RTC(j)))^N(j)-1);
    B2(2,j) = Ai(q(j))/((FC_min(j)/(PS(q(j))*RTC(q(j))))^N(q(j))-1);
end

for i = 1:N_Reles
    f(1,i) = B1(1,i) + B2(1,i);
end

Kij_geral = zeros(N_Restricoes,N_Reles);

for k = 1:4
    Kij_geral(k,s(k)) = B1(1,s(k)); % prim. max
    Kij_geral(k,u(k)) = -B1(2,s(k)); % ret. max
    Kij_geral(4+k,s(k)) = B2(1,s(k)); % prim. min
    Kij_geral(4+k,u(k)) = -B2(2,s(k)); % ret. min
end

for k = 1:5
    % Tmin
    Kij_geral(8+k,k) = -B1(1,k);
    Kij_geral(13+k,k) = -B2(1,k);
    % Tmax
    Kij_geral(18+k,k) = B1(1,k);
    Kij_geral(23+k,k) = B2(1,k);
end

%% Sistema no formato do Linprog

% F.O - Tempo mínimo de atuação dos Relés

A = Kij_geral;

[x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);

%% Tempos de Atuação dos Relês

T_pri1 = zeros(N_Reles,1);    % Pré-alocação
T_ret1 = zeros(N_Reles,1);
T_pri2 = zeros(N_Reles,1);
T_ret2 = zeros(N_Reles,1);

for i = 1:N_Reles;
    T_pri1(i,1) = x(i,1) * B1(1,i);
    T_ret1(i,1) = x(q(i),1) * B1(2,i);
    T_pri2(i,1) = x(i,1) * B2(1,i);
    T_ret2(i,1) = x(q(i),1) * B2(2,i);
end

T_coord1 = zeros(N_Reles,1);     % Configuração única!
T_coord2 = zeros(N_Reles,1);
violation1 = zeros(N_Restricoes,3);
violation2 = zeros(N_Restricoes,3);
index = 1;

for i = 1:N_Reles
    T_coord1(i,1) = p(i);
    T_coord1(i,2) = q(i);
    T_coord1(i,3) = T_pri1(i);
    T_coord1(i,4) = T_ret1(i);
    T_coord1(i,5) = T_coord1(i,4) - T_coord1(i,3);
    
    if(T_coord1(i,5) < CTI)  % Violações na Coordenação
        violation1(index,1) = p(i);
        violation1(index,2) = q(i);
        violation1(index,3) = T_coord1(i,5);
        index = index+1;
    end
end

for i = 1:N_Reles
    T_coord2(i,1) = p(i);
    T_coord2(i,2) = q(i);
    T_coord2(i,3) = T_pri2(i);
    T_coord2(i,4) = T_ret2(i);
    T_coord2(i,5) = T_coord2(i,4) - T_coord2(i,3);
    
    if(T_coord2(i,5) < CTI)  % Violações na Coordenação
        violation2(index,1) = p(i);
        violation2(index,2) = q(i);
        violation2(index,3) = T_coord2(i,5);
        index = index+1;
    end
end

t = toc;