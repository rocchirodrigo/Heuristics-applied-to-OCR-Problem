%% Sistema de 15 Barras - Tabu Search

%% Dados dos relés

load dados_15_bus.mat;
tabu_tenure_best = 1;
penalty = 1.42;
fact_FO = 10000;
fact_best_move = zeros(N_Reles,1);

%% Definições

rng('shuffle');
for i = 1:N_Reles
    PS(i,1) = PS_disp(randi(count_MAX),1);
    PS_inicial(i,vezes) = PS(i,1);
    % PS(i,1) = PS_inicial(i,vezes);        % Para PS específicos (teste)
end

options = optimoptions(@linprog,'Display','none');

%% Loop

for it = 1:N_it_MAX
    % Penalidade forte no início para soluções não factíveis, cada iteração
    % que passa há uma relaxação
        if(penalty > 1.03)
            penalty = penalty - 0.1;
        end
    for R = 1:N_Reles
        for count = 1:count_MAX
            if(tabu_short_mem(count,R) == 0)
                % Somente se o vetor não é tabu que deve-se calcular o PL!
                % Atualização do PS (incremento)
                if((R ~= 1) || (count ~= 1))
                    PS(R,1) = PS_disp(count,1);
                end
                
                % Matrizes Kij
                for j = 1:N_Reles
                    B1(1,j) = Ai/((FC_pr(j)/(PS(j)*RTC(j)))^N-1);  % Relê primário
                    B1(2,j) = Ai/((FC_bc(j)/(PS(j)*RTC(j)))^N-1);  % Relê retaguarda
                end
                
                Kij_geral = zeros(N_Restricoes,N_Reles);
                
                for k = 1:N_Restricoes
                    Kij_geral(k,p(k)) = B1(1,p(k));
                    Kij_geral(k,q(k)) = -B1(2,q(k));
                end
                
                %% Sistema no formato do Linprog
                
                % F.O - Tempo mínimo de atuação dos Relés (somente Primários)
                
                for i = 1:N_Reles
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
                
                % Exceções (bater resultado de Kida):
                A(49,13) = -Ai/((1053/(PS(23)*RTC(23)))^N-1);
                A(50,21) = -Ai/((175/(PS(21)*RTC(21)))^N-1);
                
                %% Resolução via PL
                [x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);
                PL_count = PL_count +1;
                
                if(exitflag ~= 1)
                    fval = fval * penalty;
                end
                %% Melhor movimento e atualização das Listas Tabu
                if(fval < best_fo) % Melhor movimento Local
                    best_B1 = B1;
                    best_f = f;
                    best_x = x;
                    best_move(R,1) = PS(R,1);
                    best_fo = fval;
                    pos_tabu_best = count;
                    for i = 1:N_Reles
                            position_best = PS(i,1)/0.5;
                            moves_usage(position_best,i) = moves_usage(position_best,i) + 1;
                        end
                    if(fval < best_ff)
                        % Lista dos elementos elite
                        for i = 1:N_Reles
                            position_best = PS(i,1)/0.5;
                            moves_usage(position_best,i) = moves_usage(position_best,i) + 1;
                        end
                        if(exitflag == 1)
                            fact_FO = fval;
                            fact_best_move(R,1) = PS(R,1);
                        end
                        exist_best_move(R,1) = -1; 
                        the_best_move(R,1) = PS(R,1);
                        best_ff = fval;
                    end
                else
                    if(fval < bad_fo)
                        % Caso piore a F.O mas seja a melhor opção
                        % disponível, armazenar o valor
                        bad_B1 = B1;
                        bad_f = f;
                        bad_x = x;
                        bad_move(R,1) = PS(R,1);
                        bad_fo = fval;
                    end
                end
            end
        end
        %% Melhor movimento é selecionado, listas atualizadas, valores zerados.
        % Atualização da Lista Tabu (decremeneto)
        for i = 1:count_MAX
            if(tabu_short_mem(i,R) > 0)
                tabu_short_mem(i,R) = tabu_short_mem(i,R) - 1;
            end
        end
        
        %         % Melhor movimento
        if(pos_tabu_best ~=0)
            tabu_short_mem(pos_tabu_best,R) = tabu_tenure_best;
        end
        
        % Seleção do melhor movimento
        if(exist_best_move(R,1) ~= 0)
            PS(R,1) = the_best_move(R,1);
        else PS(R,1) = best_move(R,1);
        end
        
        % Reseta todos indicadores / Contadores
        pos_tabu_best = 0;
        bad_fo = 10000;
        best_fo = 10000;
        exist_best_move(R,1) = 0;
    end
    %% Progresso das iterações
    it_results(it,1) = best_ff;
    it_results(it,2) = bad_fo_it;
    fprintf('Run = %d, it = %d, FO = %f, Tabu_Search_it = %d, Heuristic_it = %d \n',vezes, it, fact_FO, Tabu_Search_it, Heuristic_it);
    
    %% Frequência de movimentos, Aspiração / Intensificação / Diversificação
    [~,Ind1] = min(moves_usage,[],1);
    [~,Ind2] = max(moves_usage,[],1);
    
    % Intensificação - somente 1 vez senão cicla!
%     if(a ~= 0)
%         for RR = 1:N_Reles
%             if((it > 1) && (it_results(it-1,1) == it_results(it,1)))
%                 PS(RR,1) = PS_disp(Ind2(1,RR),1);
%                 moves_usage = zeros(count_MAX,N_Reles);
%             end
%         end
%     end
    
%     for RR = 1:N_Reles
%         tabu_short_mem(Ind2(1,RR),RR) = 0; % Aspiração
%     end
    clear Ind1 Ind2;
end


%% Recalcular Função Objetivo com melhores valores

% Melhores valores
PS = the_best_move;
B1 = best_B1;
f = best_f;

for k = 1:N_Restricoes
    Kij_geral(k,p(k)) = B1(1,p(k));
    Kij_geral(k,q(k)) = -B1(2,q(k));
end

%% Sistema no formato do Linprog

% F.O - Tempo mínimo de atuação dos Relés

A = Kij_geral;

% Exceções (bater resultado de Kida):
A(49,13) = -Ai/((1053/(PS(23)*RTC(23)))^N-1);
A(50,21) = -Ai/((175/(PS(21)*RTC(21)))^N-1);

[x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);
exit(vezes,1) = exitflag;


%% Tempos de Atuação dos Relês

T_pri = zeros(N_Restricoes,1);    % Pré-alocação
T_ret = zeros(N_Restricoes,1);

for i = 1:N_Reles
    T_pri(i,1) = x(i,1) * B1(1,i);
    T_ret(i,1) = x(i,1) * B1(2,i);
end

T_coord = ones(N_Restricoes,5);     % Configuração única!
violation = zeros(N_Restricoes,3);
index = 1;

for i = 1:N_Restricoes
    T_coord(i,1) = p(i);
    T_coord(i,2) = q(i);
    T_coord(i,3) = T_pri(p(i));
    T_coord(i,4) = T_ret(q(i));
    % 2 Exceções, pares R23-R13 e R24-R21
    if(i == 49 )
        T_coord(49,4) = x(q(i),1) * Ai/((1053/(PS(13)*RTC(13)))^N-1);
    end
    if(i == 50 )
        T_coord(50,4) = x(q(i),1) * Ai/((175/(PS(21)*RTC(21)))^N-1);
    end
    T_coord(i,5) = T_coord(i,4) - T_coord(i,3);
    
    if(T_coord(i,5) < CTI)  % Violações na Coordenação
        violation(index,1) = p(i);
        violation(index,2) = q(i);
        violation(index,3) = T_coord(i,5);
        index = index+1;
    end
end

t = toc;