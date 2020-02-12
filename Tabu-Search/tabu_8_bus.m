%% Sistema de 8 barras - Resolução PL

close all;
clearvars
clc
tic

%% Dados dos relés

% Plug Settings 
PS = [2.00 2.50 2.50 2.50 2.50 2.50 2.50 2.50 2.50 2.50 2.50 2.50 2.00 2.50]; % otimização Kida

f = zeros(1,14);
CTI = 0.3;
RTC = [240, 240, 160, 240, 240, 240, 160, 240, 160, 240, 240, 240, 240, 160];

% Ordem dos relês (prim-p / ret-q)
p = [1;2;2;3;4;5;6;6;7;7;8;8;9;10;11;12;12;13;14;14];

q = [6;1;7;2;3;4;5;14;5;13;7;9;10;11;12;13;14;8;1;9];

% Faltas modo normal

FC_pr = [3232, 5924, 3556, 3783, 2401, 6109, 5223, 6093, 2484, 3883, 3707, 5899, 2991, 5199];

FC_bc = [996, 3556, 2244, 2401, 1197, 3232, 1890, 2991, 1165, 2484, 2344, 3707, 987, 1874];

Icc_K1 = [FC_pr; FC_bc];

B1 = size(Icc_K1);

% Curva IEC Inversa
Ai = 0.14;
N = 0.02;

options = optimoptions(@linprog,'Display','none');

%% Matrizes Kij

for j = 1:14
    B1(1,j) = Ai/((FC_pr(j)/(PS(j)*RTC(j)))^N-1);  % Relê primário
    B1(2,j) = Ai/((FC_bc(j)/(PS(j)*RTC(j)))^N-1);  % Relê retaguarda
end

Kij_geral = zeros(20,14);

for k = 1:20
        Kij_geral(k,p(k)) = B1(1,p(k));
        Kij_geral(k,q(k)) = -B1(2,q(k));
end

%% Sistema no formato do Linprog

% F.O - Tempo mínimo de atuação dos Relés (somente Primários)

 for i = 1:14
     f(1,i) = B1(1,i);
 end
 
% Valores Mín/Máx de Dial dos Relés R1-R6
lb = ones(1,14);
lb = lb * 0.1;
ub = ones(1,14);
ub = ub * 1.1;

% Lado direito da Desigualdade
b = ones(20,1);
b = b * -CTI;

A = Kij_geral;

% Exceções (sem exceções neste sistema):

%% Resolução via PL

[x,fval] = linprog(f,A,b,[],[],lb,ub,[],options);

%% Tempos de Atuação dos Relês

T_pri = zeros(20,1);    % Pré-alocação
T_ret = zeros(20,1);

for i = 1:14
    T_pri(i,1) = x(i,1) * B1(1,i);
    T_ret(i,1) = x(i,1) * B1(2,i);
end

T_coord = ones(20,5);     % Configuração única!

for i = 1:20
    T_coord(i,1) = p(i);
    T_coord(i,2) = q(i);
    T_coord(i,3) = T_pri(p(i));
    T_coord(i,4) = T_ret(q(i));
    T_coord(i,5) = T_ret(q(i)) - T_pri(p(i));
end

t = toc;