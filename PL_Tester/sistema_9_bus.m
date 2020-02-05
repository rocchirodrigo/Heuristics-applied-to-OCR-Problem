%% Sistema de 9 barras - Resolução PL

close all;
clearvars
clc
tic

%% Dados dos relés

% Plug Settings 
PS = [1.81 0.69 1.49 1.39 0.94 1.63 1.63 0.94 1.39 0.90 0.87 0.81 1.37 1.39 1.55 0.99 1.36 1.55 1.12 0.76 1.60 0.21 1.43 1.38]'; % otimização Kida

f = zeros(1,24);
CTI = 0.2;      % Intervalo de Coordenação
%Tmin_pr = 0.2;    % Tempo mínimo de atuação dos relês primários
RTC = 500;

% Ordem dos relês (prim-p / ret-q)
p = [1;1;2;3;4;5;6;6;7;7;8;9;10;11;12;12;13;13;14;14;15;15;16;16;18;18;20;20;22;22;24;24];

q = [15;17;4;1;6;3;8;23;5;23;10;7;12;9;14;21;11;21;16;19;13;19;2;17;2;15;13;16;11;14;5;8];

% Faltas modo normal
% Obs: utilizado Icc(13) = 1503 [A] para o par R12-R13 e 1053 [A] para o
% par R23-R13! No par R24-R21 foi utilizado 1326 [A]!

FC_pr = [4863.6, 1634.4, 2811.4, 2610.5, 1778.0, 4378.5, 4378.5, 1778.0, 2610.5, 2811.4, 1634.4, 2811.4, 3684.5, 4172.5, 4172.5, 3684.5, 7611.2, ...
    2271.7, 7435.8, 2624.2, 7611.2, 2271.7, 7914.7, 1665.5];

FC_bc = [1361.6, 653.6, 1124.4, 1044.2, 711.2, 1226.0 1226.0, 711.2, 1044.2 1124.4 653.6, 787.2, 1031.7, 1168.3, 1168.3, 1031.7, 1293.9, 1953.7, ...
     1264.1, 2256.8, 1293.9, 1953.7, 1345.5, 1432.3];

Icc_K1 = [FC_pr; FC_bc];

B1 = size(Icc_K1);

% Curva IEC Inversa
Ai = 0.14;
N = 0.02;

options = optimoptions(@linprog,'Display','none');

%% Matrizes Kij

for j = 1:24
    B1(1,j) = Ai/((FC_pr(j)/(PS(j)*RTC))^N-1);  % Relê primário
    B1(2,j) = Ai/((FC_bc(j)/(PS(j)*RTC))^N-1);  % Relê retaguarda
end

Kij_geral = zeros(56,24);

for k = 1:56
    if(k <=32)  % Restrições de Coordenação (par primário-retaguarda)
        Kij_geral(k,p(k)) = B1(1,p(k));
        Kij_geral(k,q(k)) = -B1(2,q(k));    
    else Kij_geral(k,k-32) = -B1(1,k-32);   % Restrição Tmin = 0.2s
    end
end

%% Sistema no formato do Linprog

% F.O - Tempo mínimo de atuação dos Relés (somente Primários)

 for i = 1:24
     f(1,i) = B1(1,i);
 end
 
% Valores Mín/Máx de Dial dos Relés R1-R6
lb = ones(1,24) * 0.025;
ub = ones(1,24) * 1.2;

% Lado direito da Desigualdade
b = ones(56,1) * -CTI;

A = Kij_geral;

% Exceções (sem exceções neste sistema):

%% Resolução via PL

[x,fval] = linprog(f,A,b,[],[],lb,ub,[],options);

%% Tempos de Atuação dos Relês

T_pri = zeros(24,1);    % Pré-alocação
T_ret = zeros(24,1);

for i = 1:24
    T_pri(i,1) = x(i,1) * B1(1,i);
    T_ret(i,1) = x(i,1) * B1(2,i);
end

T_coord = ones(32,5);     % Configuração única!

for i = 1:32
    T_coord(i,1) = p(i);
    T_coord(i,2) = q(i);
    T_coord(i,3) = T_pri(p(i));
    T_coord(i,4) = T_ret(q(i));
    T_coord(i,5) = T_ret(q(i)) - T_pri(p(i));
end

t = toc;