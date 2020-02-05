%% Sistema de 6 barras - Resolução PL

close all;
clearvars
clc
tic

%% Dados dos relés

% Plug Settings 
PS = [1.5000 1.5000 1.2776 1.5000 1.2500 1.3808 1.2500 1.2500 1.2500 1.3968 1.5000 1.5000 1.4329 1.500]'; % otimização Kida

f = zeros(1,14);
CTI = 0.2;
% Não tem RTC!! = [0.2585,0.2585,0.4863,0.4863,0.7138,0.7138,1.746,1.746,1.0424,1.0424,0.7729,0.7729,0.5879,0.5879];

% Ordem dos relês para as restrições
p = [1 1 1 2 2 3 3 3 4 4 5 5 5 6 6 7 7 7 7 9 9 9 9 11 11 11 11 12 12 12 12 13 13 13 14 14 14 14]; % Kida
%p=[1 1 2 3 3 3 4 5 5 6 7 7 9 9 11 11 12 12 13 13 14 14 1 2 4 5 6 7 7 9 9 11 11 12 12 13 14 14];
q = [8 8 11 3 3 10 10 13 1 1 13 13 14 3 3 2 2 11 11 4 4 13 13 6 6 14 14 2 2 8 8 6 12 12 4 4 10 10]; % Kida
%q=[8 11 3 10 10 13 1 12 14 3 2 11 4 13 6 14 2 8 6 12 4 10 8 3 1 12 3 2 11 4 13 6 14 2 8 12 4 10];

% Ordem dos tempos de atuação no início / fim da linha
r =[1 2 3 4 5 6 7 8 9 10 11 12 13 14];

s = [2 1 4 3 6 5 8 7 10 9 12 11 14 13];

% Faltas (somente FC_pr e FC_bc são utilizadas para restrições!)

%FC_pr = [20.79 20.79 10.59 6.11 9.44 9.44 8.53 5.96 5.96 8.73 2.39 2.39 4.46 4.46 3.34 3.34 7.91 7.91 7.37 4.11 9.11 9.11 9.79 23.02 13.7 2.74 8.73 2.20 2.20 5.06 5.06 5.05 5.05 4.92 4.92 7.37 4.93 4.93];   
FC_pr = [20.794 20.794 10.590 6.112 9.438 9.438 8.529 5.964 5.96 8.73 2.39 2.39 4.46 4.46 3.34 3.34 7.91 7.91 7.37 4.11 9.11 9.11 9.79 23.02 13.7 2.74 8.73 2.20 2.20 5.06 5.06 5.05 5.05 4.92 4.92 7.37 4.93 4.93];   

FC_bc = [2.34 1.67 1.28 2.46 0.09 2.55 3.43 3.29 2.92 3.01 7.87 2.77 6.24 2.75 1.55 2.51 5.97 2.62 2.25 2.38 5.31 2.67 1.68 3.43 5.90 1.88 2.31 7.24 2.55 7.07 3.12 2.54 3.55 1.83 1.91 3.53 1.80 1.94];

FC_ini = [9.7915 10.5903 6.1121 8.5291 2.7382 3.8776 2.2006 3.2176 4.4645 3.3827 3.3433 4.9173 4.1067 9.1072];

FC_fim = [23.0155 20.7938 13.7037 9.4380 8.7342 5.9643 3.6480 2.3931 3.7126 5.0553 7.9110 5.0466 4.9347 7.3737];

Icc_K1 = [FC_pr; FC_bc];
Icc_K2 = [FC_ini; FC_fim];

B1 = size(Icc_K1);
B2 = size(Icc_K2);

% Curva IEC Inversa
Ai = 0.14;
N = 0.02;

options = optimoptions(@linprog,'Display','none');

%% Matrizes Kij

% Restrições de Coordenação
for j = 1:14
    B1(1,j) = Ai/((FC_pr(j)/(PS(j)))^N-1);  % Relê primário
    B1(2,j) = Ai/((FC_bc(j)/(PS(j)))^N-1);  % Relê retaguarda
end

% Restrições de Tmín/Tmáx
for j = 1:14
    B2(1,j) = Ai/((FC_ini(j)/(PS(j)))^N-1);  % Relê início da linha
    B2(2,j) = Ai/((FC_fim(j)/(PS(j)))^N-1);  % Relê fim da linha
end

Kij_geral = zeros(94,14);

for k = 1:38    % Restrições Tret-Tprim
        Kij_geral(k,p(k)) = B1(1,p(k));
        Kij_geral(k,q(k)) = -B1(2,q(k));
end

for k = 39:52   % Restrições Tmín
    Kij_geral(k,k-38) = -B2(1,r(k-38));
end

for k = 53:66   % Restrições Tmáx
    Kij_geral(k,k-52) = B2(1,r(k-52));
end

for k = 67:80   % Restrições Tmín
    Kij_geral(k,k-66) = -B2(2,k-66);
end

for k = 81:94   % Restrições Tmáx
    Kij_geral(k,k-80) = B2(2,k-80);
end


%% Sistema no formato do Linprog

% F.O - Tempo mínimo de atuação dos Relés (somente Primários)

 for i = 1:14
     f(1,i) = B1(1,i);
 end
 
% Valores Mín/Máx de Dial dos Relés R1-R6
lb = ones(1,14) * 0.05;
ub = ones(1,14) * 1.1;

% Lado direito da Desigualdade
b = ones(94,1);

for i = 1:38
    b(i) = -CTI;
end

for i = 39:52
    b(i) = -0.05;
end

for i = 53:66
    b(i) = 1.00;
end

for i = 67:80
    b(i) = -0.05;
end

for i = 81:94
    b(i) = 1.00;
end

A = Kij_geral;

% Exceções (bater resultado de Kida):
%A(5,p(5)) = 0;
%A(5,q(5)) = 0;
%b(5) = 1000;

%% Resolução via PL

[x,fval] = linprog(f,A,b,[],[],lb,ub,[],options);

%% Tempos de Atuação dos Relês

T_pri = zeros(38,1);    % Pré-alocação
T_ret = zeros(38,1);
T_ini = zeros(38,1);
T_fim = zeros(38,1);

for i = 1:14
    T_pri(i,1) = x(i,1) * B1(1,i);
    T_ret(i,1) = x(i,1) * B1(2,i);
    T_ini(i,1) = x(i,1) * B2(1,i);
    T_fim(i,1) = x(i,1) * B2(2,i);
end

T_coord = ones(38,7);     % Configuração única!

for i = 1:38
    T_coord(i,1) = p(i);
    T_coord(i,2) = q(i);
    T_coord(i,3) = T_pri(p(i));
    T_coord(i,4) = T_ret(q(i));
    T_coord(i,5) = T_ret(q(i)) - T_pri(p(i));
    if(i<=14)
    T_coord(i,6) = T_ini(i);
    T_coord(i,7) = T_fim(s(i));
    else T_coord(i,6) = 0; T_coord(i,7) = 0;
    end
end

t = toc;