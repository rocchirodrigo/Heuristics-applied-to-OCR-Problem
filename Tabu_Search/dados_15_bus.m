%% Arquivo de Dados para Sistema de 15 Barras

N_Reles = 42;   N_Restricoes = 82; count_MAX = 5;

% Plug Settings
PS_disp = zeros(count_MAX,1);

for i = 1:count_MAX
    PS_disp(i,1) = 0.5 + 0.5*(i-1);
end

% Geração de Ponto Inicial Aleatório
PS = zeros(N_Reles,1);

% Flags
PL_count = 0;

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
f = zeros(1,N_Reles);

% Curva IEC Inversa
Ai = 0.14;
N = 0.02;

%% Vetores Tabu

% Critérios de Frequência
moves_usage = zeros(count_MAX, N_Reles);      % Valores de PS mais utilizados (Diversificação)
tabu_short_mem = zeros(count_MAX, N_Reles);           % valores de PS proibidos, dimensão (PS_disp x Relês)

% Melhor movimento
exist_flag = 0;                         % Existe melhor movimento disponível?
best_move = ones(N_Reles,1) * 0.5;
the_best_move = ones(N_Reles,1) * 0.5;
the_best_move = zeros(N_Reles,1);
best_B1 = size(B1);
best_fo = 10000;                        % melhor resultado
best_ff = 10000;
best_fo_it = 10000;                     % para definir parada
pos_tabu_best = 0;                      % posição do melhor movimento
position_best = 0;

% Bad Move (Diversificação apenas!)
exist_best_move = zeros(N_Reles,1);          % Indica se há um melhor movimento disponível
bad_move = ones(N_Reles,1) * 0.5;            % Movimento que menos piora a F.O
bad_fo = 1000;
bad_fo_it = 0;
bad_B1 = size(B1);

% Iterações
%it_max = 20;
it_results = zeros(200,2);              % F.O
it_lim_div = 0;                         % Iteração máx. estágio de Div.
it_lim_it = 0;                          % Iteração máx. estágio de Int.
Tabu_Search_it = 0;                     % Iteração na qual foi pra Diversificação
Heuristic_it = 1;
Intensification_it = 0;

% Duração do banimento
tabu_tenure_best = 2;                   % banimento do melhor movimento
penalty = 1.10;