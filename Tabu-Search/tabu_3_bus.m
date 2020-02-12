%% Sistema de 3 barras - Tabu Search

close all;
clearvars
clc
tic

%% Dados dos rel�s

% Plug Settings 
%PS = [2.5 2.0 3.0 2.5 2.5 1.5]; % Otimiza��o Kida
PS_disp = zeros(8,1);
PS = ones(6,1) * 0.5;

for i = 1:8
        PS_disp(i,1) = 1.5 + 0.5*(i-1);
end

% Flags
PL_count = 0;
flag = 0;

CTI = 0.2;
RTC = [300/5, 200/5, 200/5, 300/5, 200/5, 400/5];

% Ordem dos rel�s (prim-p / ret-q)
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
f = zeros(1,6); 
best_f = zeros(1,6);

% Curva IEC Inversa
Ai = 0.14;
N = 0.02;

options = optimoptions(@linprog,'Display','none');

%% Vetores Tabu

tabu_short_mem = zeros(8,6);    % valores de PS proibidos, dimens�o (PS_disp x Rel�s)
tabu_best_sol = zeros();        % melhores valores de PS para diversifica��o
best_move = zeros(6,1);         % melhor movimento (no caso PS que n�o seja tabu)
best_fo = 10000;                % melhor resultado
best_fo_it = 10000;             % para definir parada
it_max = 4;
it_results = zeros(20,1);      % F.O
pos_tabu = 0;                   % posi��o do melhor movimento
tabu_state = 0;                 % Est�gio do algoritmo:
                                %   0: Intensifica��o (busca)
                                %   1: Diversifica��o
                                %   2: Resultado final
for it = 1:it_max
    switch tabu_state
        case 0
%% Intensifica��o
            for R = 1:6
                for count = 1:8
                    if(tabu_short_mem(count,R) == 0)    % Somente se o vetor n�o � tabu
                        % que deve-se calcular o PL!
                        
                        % Atualiza��o do PS (incremento)
                        PS(R,1) = PS_disp(count,1);
                        
                        %% Matrizes Kij
                        for j = 1:6
                            % Matriz faltas normais
                            B1(1,j) = Ai/((FC_pr(j)/(PS(j)*RTC(j)))^N-1);
                            B1(2,q(j)) = Ai/((FC_bc(j)/(PS(q(j))*RTC(q(j))))^N-1);
                            % Matriz transiente
                            B2(1,j) = Ai/((FC_pr_Tr(j)/(PS(j)*RTC(j)))^N-1);
                            B2(2,q(j)) = Ai/((FC_bc_Tr(j)/(PS(q(j))*RTC(q(j))))^N-1);
                        end
                        
                        Kij_geral = zeros(12,6);
                        
                        for k = 1:6
                            Kij_geral(k,p(k)) = B1(1,p(k));
                            Kij_geral(k,q(k)) = -B1(2,q(k));
                            Kij_geral(k+6,p(k)) = B2(1,p(k));
                            Kij_geral(k+6,q(k)) = -B2(2,q(k));
                        end
                        
                        %% Sistema no formato do Linprog
                        
                        % F.O - Tempo m�nimo de atua��o dos Rel�s
                        
                        for i = 1:6
                            f(1,i) = B1(1,i);
                        end
                        
                        % Valores M�n/M�x de Dial dos Rel�s R1-R6
                        lb = ones(1,6) * 0.1;
                        ub = ones(1,6) * 1.1;
                        
                        % Lado direito da Desigualdade
                        b = ones(12,1) * -CTI;
                        
                        A = Kij_geral;
                        
                        %% Resolu��o via PL
                        
                        [x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);
                        PL_count = PL_count + 1;
                        
                        if(fval < best_fo)
                            best_f = f;
                            best_x = x;
                            best_move(R,1) = PS(R,1);
                            best_fo = fval;
                            pos_tabu = count;
                        end
                        
                        % Seleciona o melhor movimento e o pro�be de ser
                        % executado novamente
                        if(count == 8)      
                            PS(R,1) = best_move(R,1);
                            tabu_short_mem(pos_tabu,R) = 4;
                        end
                    end
                end
            end
            
            % Atualiza��o da Lista Tabu (decremento)
            for i = 1:count
                for j = 1:R
                    if(tabu_short_mem(i,j) > 0)
                        tabu_short_mem(i,j) = tabu_short_mem(i,j) - 1;
                    end
                end
            end
            
            % H� melhoras na Fun��o Objetivo?
            if(best_fo < best_fo_it)
                best_fo_it = best_fo;
            else tabu_state = 2;
            end

%% Diversifica��o            
        case 1
            flag = 1;
            
%% Recalcular Fun��o Objetivo com melhores valores            
        case 2            
            flag = 2;
            PS = best_move;
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
                        
                        % F.O - Tempo m�nimo de atua��o dos Rel�s
                        
                        for i = 1:6
                            f(1,i) = B1(1,i);
                        end

                        A = Kij_geral;
                        
                        [x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);
                        break;
    end
    a = 1;
    
    % Progresso das itera��es
    it_results(it,1) = best_fo;
end

%% Tempos de Atua��o dos Rel�s

T_pri1 = zeros(1,6);    % Pr�-aloca��o
T_ret1 = zeros(1,6);
T_pri2 = zeros(1,6); 
T_ret2 = zeros(1,6);

for i = 1:6
    T_pri1(1,i) = x(i,1) * B1(1,i);
    T_ret1(1,i) = x(i,1) * B1(2,i);
    T_pri2(1,i) = x(i,1) * B2(1,i);
    T_ret2(1,i) = x(i,1) * B2(2,i);
    
end

T_coord1 = ones(6,1);     % Configura��o �nica!
T_coord2 = ones(6,1);

for i = 1:6
    T_coord1(i,1) = T_ret1(q(i)) - T_pri1(p(i));
    T_coord2(i,1) = T_ret2(q(i)) - T_pri2(p(i));
end

t = toc;