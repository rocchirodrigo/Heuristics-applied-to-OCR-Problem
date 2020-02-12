%% Sistema de 3 barras - Resolução PL

close all;
% clearvars
clc
tic

%% Dados dos relés

% Plug Settings 
%PS = [1.50 1.00 2.00 1.00 2.00 2.00 2.00 1.50 2.00 1.50 1.50 1.50 2.00 1.00 1.00 1.50 2.00 1.00 2.00 1.50 0.50 1.50 1.00 1.50 2.00 1.50 2.00 2.50 1.50 2.00 2.00 1.50 2.50 2.50 2.00 2.00 2.50 2.50 2.50 2.50 2.50 1.50]; % otimização Kida
% PS = ones(1,42) * 0.5;

N_Reles = 42;
N_Restricoes = 82;

% PS = ones(N_Reles,1) * 0.5; 
f = zeros(1,42);
CTI = 0.2;
RTC=[160, 240, 160, 240, 160, 120, 120, 240, 120, 160, 240, 240, 160, 240, 240, 120, 80, 320, 160, 320, 320, ...
    80, 240, 120, 120, 120, 120, 120, 320, 80, 120, 120, 120, 80, 120, 160, 160, 80, 80, 160, 80, 160];

% Ordem dos relês (prim-p / ret-q)
p=[1;2;2;3;3;4;4;4;5;6;6;7;7;8;8;8;9;9;10;11;11;11;12;12;13;14;14;15;15;16;16;17;17;18;18;18;19;19;19;20;20;20;21;21;21;
    22;22;23;23;24;24;25;25;26;26;27;27;28;28;29;29;29;30;30;31;31;32;32;33;33;34;34;35;35;36;37;38;39;40;41;41;42
    ];

q=[6;4;16;1;16;7;12;20;2;8;10;5;10;3;12;20;5;8;14;3;7;20;13;24;9;11;24;1;4;18;26;15;26;19;22;30;3;7;12;17;22;30;17;19;30;
    23;34;11;13;21;34;15;18;28;36;25;36;29;32;17;19;22;27;32;27;29;33;42;21;23;31;42;25;28;38;35;40;37;41;31;33;39
    ];

% Faltas modo normal
% Obs: utilizado Icc(13) = 1503 [A] para o par R12-R13 e 1053 [A] para o
% par R23-R13! No par R24-R21 foi utilizado 1326 [A]!

FC_pr = [3621, 4597, 3984, 4382, 3319, 2647, 2497, 4695, 2943, 3568, 4342, 4195, 3402, 4606, 4712, 2225, 1875, 8426, 3998, 7662, 8384, ...
    1950, 4910, 2296, 2289, 2300, 2011, 2525, 8346, 1736, 2867, 2069, 2305, 1715, 2095, 3283, 3301, 1403, 1434, 3140, 1971, 3295]; 

FC_bc = [853, 922, 1424, 1477, 1397, 1233, 1111, 1548, 1009, 1100, 1475, 1463, 1503, 1175, 969, 743, 599, 1320, 1372, 1808, 1326, ...
    642, 979, 753, 903, 905, 1039, 1192, 1828, 681, 809, 697, 1162, 970, 910, 1109, 1434, 882, 896, 1403, 745, 907];

Icc_K1 = [FC_pr; FC_bc];

B1 = size(Icc_K1);

% Curva IEC Inversa
Ai = 0.14;
N = 0.02;

options = optimoptions(@linprog,'Display','none');

%% Matrizes Kij

for j = 1:42
    B1(1,j) = Ai/((FC_pr(j)/(PS(j)*RTC(j)))^N-1);  % Relê primário
    B1(2,j) = Ai/((FC_bc(j)/(PS(j)*RTC(j)))^N-1);  % Relê retaguarda
end

Kij_geral = zeros(82,42);

for k = 1:82
        Kij_geral(k,p(k)) = B1(1,p(k));
        Kij_geral(k,q(k)) = -B1(2,q(k));
end

%% Sistema no formato do Linprog

% F.O - Tempo mínimo de atuação dos Relés (somente Primários)

 for i = 1:42
     f(1,i) = B1(1,i);
 end
 
% Valores Mín/Máx de Dial dos Relés R1-R6
lb = ones(1,42);
lb = lb * 0.1;
ub = ones(1,42);
ub = ub * 1.1;

% Lado direito da Desigualdade
b = ones(82,1);
b = b * -CTI;

A = Kij_geral;

% Exceções (bater resultado de Kida):
A(49,13) = -Ai/((1053/(PS(23)*RTC(23)))^N-1);
A(50,21) = -Ai/((175/(PS(21)*RTC(21)))^N-1);
%% Resolução via PL

[x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);

%% Tempos de Atuação dos Relês

T_pri = zeros(82,1);    % Pré-alocação
T_ret = zeros(82,1);

for i = 1:42
    T_pri(i,1) = x(i,1) * B1(1,i);
    T_ret(i,1) = x(i,1) * B1(2,i);
end

T_coord = ones(82,5);     % Configuração única!

for i = 1:82
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
end

t = toc;