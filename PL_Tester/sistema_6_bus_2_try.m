%% Sistema de 6 barras - Resolu��o PL

close all;
clearvars
clc
tic

%% Dados dos rel�s

% Plug Settings 
PS = [1.5000 1.5000 1.2776 1.5000 1.2500 1.3808 1.2500 1.2500 1.2500 1.3968 1.5000 1.5000 1.4329 1.500]'; % otimiza��o Kida

f = zeros(1,14);
CTI = 0.2;
% N�o tem RTC!

% Ordem dos rel�s para as restri��es
p = [1 1 1 2 2 3 3 3 4 4 5 5 5 6 6 7 7 7 7 9 9 9 9 11 11 11 11 12 12 12 12 13 13 13 14 14 14 14]; % Kida

q = [8 8 11 3 3 10 10 13 1 1 13 13 14 3 3 2 2 11 11 4 4 13 13 6 6 14 14 2 2 8 8 6 12 12 4 4 10 10]; % Kida

% Ordem dos tempos de atua��o no in�cio / fim da linha
r =[1 2 3 4 5 6 7 8 9 10 11 12 13 14];

s = [2 1 4 3 6 5 8 7 10 9 12 11 14 13];

% Faltas (somente FC_pr e FC_bc s�o utilizadas para restri��es!)

FC_pr = [20.794 9.791 20.794 10.590 23.015 6.112 9.438 9.438 8.529 13.704 5.964 2.738 5.964 8.734 8.734 2.393 2.201 2.393 2.201 4.465 5.055 4.465 5.055 3.343 5.047 3.343 5.047 7.911 4.917 7.911 4.917 7.374 4.107 7.374 9.107 4.935 9.107 4.935];
FC_bc = [];

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

% Restri��es de Coordena��o
for j = 1:14
    B1(1,j) = Ai/((FC_pr(j)/(PS(j)))^N-1);  % Rel� prim�rio
    B1(2,j) = Ai/((FC_bc(j)/(PS(j)))^N-1);  % Rel� retaguarda
end

% Restri��es de Tm�n/Tm�x
for j = 1:14
    B2(1,j) = Ai/((FC_ini(j)/(PS(j)))^N-1);  % Rel� in�cio da linha
    B2(2,j) = Ai/((FC_fim(j)/(PS(j)))^N-1);  % Rel� fim da linha
end

Kij_geral = zeros(94,14);

for k = 1:38    % Restri��es Tret-Tprim
        Kij_geral(k,p(k)) = B1(1,p(k));
        Kij_geral(k,q(k)) = -B1(2,q(k));
end

for k = 39:52   % Restri��es Tm�n
    Kij_geral(k,k-38) = -B2(1,r(k-38));
end

for k = 53:66   % Restri��es Tm�x
    Kij_geral(k,k-52) = B2(1,r(k-52));
end

for k = 67:80   % Restri��es Tm�n
    Kij_geral(k,k-66) = -B2(2,k-66);
end

for k = 81:94   % Restri��es Tm�x
    Kij_geral(k,k-80) = B2(2,k-80);
end


%% Sistema no formato do Linprog

% F.O - Tempo m�nimo de atua��o dos Rel�s (somente Prim�rios)

 for i = 1:14
     f(1,i) = B1(1,i);
 end
 
% Valores M�n/M�x de Dial dos Rel�s R1-R6
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

% Exce��es (bater resultado de Kida):
%A(5,p(5)) = 0;
%A(5,q(5)) = 0;
%b(5) = 1000;

%% Resolu��o via PL

[x,fval] = linprog(f,A,b,[],[],lb,ub,[],options);

%% Tempos de Atua��o dos Rel�s

T_pri = zeros(38,1);    % Pr�-aloca��o
T_ret = zeros(38,1);
T_ini = zeros(38,1);
T_fim = zeros(38,1);

for i = 1:14
    T_pri(i,1) = x(i,1) * B1(1,i);
    T_ret(i,1) = x(i,1) * B1(2,i);
    T_ini(i,1) = x(i,1) * B2(1,i);
    T_fim(i,1) = x(i,1) * B2(2,i);
end

T_coord = ones(38,7);     % Configura��o �nica!

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