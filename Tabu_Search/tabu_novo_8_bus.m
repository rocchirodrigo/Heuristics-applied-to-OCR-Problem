%% Sistema de 15 Barras - Tabu Search

%% Dados dos relés

load dados_8_bus.mat;

%% Definições

rng('shuffle');
for i = 1:N_Reles
    PS(i,1) = PS_disp(randi(count_MAX),1);
    PS_inicial(i,vezes) = PS(i,1);
%     PS(i,1) = PS_inicial(i,vezes);        % Para PS específicos (teste)
end

options = optimoptions(@linprog,'Display','none');

%% Loop

for it = 1:N_it_MAX
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
                
                %% Resolução via PL
                [x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);
                PL_count = PL_count +1;
                
                %% Melhor movimento e atualização das Listas Tabu
                if(exitflag == 1)
                    if(fval <= best_fo) % Melhor movimento Local
                        best_B1 = B1;
                        best_f = f;
                        best_x = x;
                        best_move(R,1) = PS(R,1);
                        best_fo = fval;
                        if(fval <= best_ff) % Seleção do melhor movimento geral
                            best_ff = fval;
                        end
                        pos_tabu_best = count;
                        exist_best_move(R,1) = -1;
                        moves_usage(count,R) = moves_usage(count,R) +1;
                    else
                        if(fval <= bad_fo)
                            % Caso piore a F.O mas seja a melhor opção
                            % disponível, armazenar o valor
                            bad_B1 = B1;
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
                    
                    % Verificar se precisa disso!!!!!
                    PS(R,1) = best_move(R,1);
                    
                    % Melhor movimento
                    if(pos_tabu_best ~=0)
                        tabu_short_mem(pos_tabu_best,R) = tabu_tenure_best;
                    end
                    
                    % 1ª Iteração: Heurística dá solução inicial
                    % 2ª Iteração: Busca Tabu com banimentos
                    if(it > 1)
                        best_fo = 10000;
                        Heuristic_it = 0;
                        if(Intensification_it == 0)
                            tabu_tenure_best = 2;
                            Tabu_Search_it = 1;
                        end
                    else tabu_tenure_best = 0;
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
    % Progresso das iterações
    it_results(it,1) = best_ff;
    it_results(it,2) = bad_fo_it;
    fprintf('Run = %d, it = %d, FO = %f, Tabu_Search_it = %d, Heuristic_it = %d \n',vezes, it, best_ff, Tabu_Search_it, Heuristic_it);
    
    % Algoritmo estagnado, faz Intensificação
    if((it > Intensification_it + 2) && (Intensification_it == 0))
        if((it > 2) && (it_results(it-1,1) == it_results(it,1)))
            Intensification_it = it;
            Tabu_Search_it = 2;
            % Selecionar movimento mais frequente.
            % Na intensificação, não se usam banimentos
            [M,Ind] = min(moves_usage,[],1);
            for R = 1:N_Reles
                PS(R,1) = PS_disp(Ind(1,R),1);
            end
            clear M Ind;
%             moves_usage = zeros(count_MAX, N_Reles);
            tabu_short_mem = zeros(count_MAX, N_Reles);
            tabu_tenure_best = 0;
        end
    end
end


%% Recalcular Função Objetivo com melhores valores

% Melhores valores
PS = best_move;
B1 = best_B1;
f = best_f;

for k = 1:N_Restricoes
    Kij_geral(k,p(k)) = B1(1,p(k));
    Kij_geral(k,q(k)) = -B1(2,q(k));
end

%% Sistema no formato do Linprog

% F.O - Tempo mínimo de atuação dos Relés

A = Kij_geral;

[x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);



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