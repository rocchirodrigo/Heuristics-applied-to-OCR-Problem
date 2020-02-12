%% Sistema de 3 barras - Resolução PL

close all;
clearvars
clc
tic

%% Dados dos relés

% Plug Settings 
% PS = [2.5 2.0 3.0 2.5 2.5 1.5]; % Otimização Kida
PS = [1.5 1.5 1.5 1.5 1.5 1.5]; % Ponto inicial
step = 0.5;

f = zeros(1,6);
f_best = 10000;
fo_best = zeros(1,6);
x_best = zeros(6,1);
PS_best = zeros(1,6);

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

% Curva IEC Inversa
Ai = 0.14;
N = 0.02;

options = optimoptions(@linprog,'Display','none');

%% Loop

for R = 1:6
    for count = 1:8

% Matrizes Kij

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

% F.O - Tempo mínimo de atuação dos Relés

 for i = 1:6
     f(1,i) = B1(1,i);
 end
 
% Valores Mín/Máx de Dial dos Relés R1-R6
lb = ones(1,6) * 0.1;
ub = ones(1,6) * 1.1;

% Lado direito da Desigualdade
b = ones(12,1) * -CTI;

A = Kij_geral;

%% Resolução do PL Relaxado

[x,fval] = linprog(f,A,b,[],[],lb,ub,[],options);
PL_count = PL_count +1;

%% Algoritmo Colosso
        if(fval <= f_best)  % Armazena a melhor solução e incrementa
            x_best = x;
            f_best = fval;
            fo_best = f;
            PS_best = PS;
            if(count < 8)
                PS(1,R) = PS(1,R) + step;
            end
        else % Decrementa caso não seja a melhor solução
                PS(1,R) = PS(1,R) - step;

            if(R == 6)
                x = x_best;
                f = fo_best;
                fval = f_best;
                PS = PS_best;
            else PS(1,R+1) = PS(1,R+1) + step; 
            end
             break;
        end
    end
end

%% Tempos de Atuação dos Relês

T_pri1 = zeros(1,6);    % Pré-alocação
T_ret1 = zeros(1,6);
T_pri2 = zeros(1,6); 
T_ret2 = zeros(1,6);

for i = 1:6
    T_pri1(1,i) = x(i,1) * B1(1,i);
    T_ret1(1,i) = x(i,1) * B1(2,i);
    T_pri2(1,i) = x(i,1) * B2(1,i);
    T_ret2(1,i) = x(i,1) * B2(2,i);
    
end

T_coord1 = ones(6,3);     % Configuração única!
T_coord2 = ones(6,3);

for i = 1:6
    T_coord1(i,1) = T_pri1(p(i));
    T_coord1(i,2) = T_ret1(q(i));
    T_coord1(i,3) = T_ret1(q(i)) - T_pri1(p(i));
    T_coord2(i,1) = T_pri2(p(i));
    T_coord2(i,2) = T_ret2(q(i));
    T_coord2(i,3) = T_ret2(q(i)) - T_pri2(p(i));
end

t = toc;