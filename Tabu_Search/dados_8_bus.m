%% Sistema de 8 Barras - Tabu Search

%% Dados dos relés
N_Reles = 14;   N_Restricoes = 20; count_MAX = 7;

% Plug Settings

PS_disp = [0.5 0.6 0.8 1.0 1.5 2.0 2.5]';

% Geração de Ponto Inicial Aleatório
PS = zeros(N_Reles,1);

% Flags
PL_count = 0;

CTI = 0.3;
RTC = [240, 240, 160, 240, 240, 240, 160, 240, 160, 240, 240, 240, 240, 160];

% Ordem dos relês (prim-p / ret-q)
p = [1;2;2;3;4;5;6;6;7;7;8;8;9;10;11;12;12;13;14;14];

q = [6;1;7;2;3;4;5;14;5;13;7;9;10;11;12;13;14;8;1;9];

% Níveis de Corrente de Falta.
% Sem exceções!

FC_pr = [3232, 5924, 3556, 3783, 2401, 6109, 5223, 6093, 2484, 3883, 3707, 5899, 2991, 5199];

FC_bc = [996, 3556, 2244, 2401, 1197, 3232, 1890, 2991, 1165, 2484, 2344, 3707, 987, 1874];

Icc_K1 = [FC_pr; FC_bc];

B1 = size(Icc_K1);
f = zeros(1,N_Reles);

% Curva IEC Inversa
Ai = 0.14;
N = 0.02;

options = optimoptions(@linprog,'Display','none');

%% Vetores Tabu

% Critérios de Frequência
moves_usage = zeros(count_MAX, N_Reles);      % Valores de PS mais utilizados (Diversificação)
tabu_short_mem = zeros(count_MAX, N_Reles);           % valores de PS proibidos, dimensão (PS_disp x Relês)

% Melhor movimento
exist_flag = 0;                         % Existe melhor movimento disponível?
best_move = ones(N_Reles,1) * 0.5;
best_B1 = size(B1);
best_fo = 10000;                        % melhor resultado
best_ff = 10000;
best_fo_it = 10000;                     % para definir parada
pos_tabu_best = 0;                      % posição do melhor movimento

% Bad Move (Diversificação apenas!)
% exist_best_move = zeros(N_Reles,1);          % Indica se há um melhor movimento disponível
bad_move = ones(N_Reles,1) * 0.5;            % Movimento que menos piora a F.O
bad_fo = 1000;
bad_fo_it = 0;
bad_B1 = size(B1);

% Iterações
%it_max = 20;
it_results = zeros(200,2);              % F.O
it_lim_div = 0;                         % Iteração máx. estágio de Div.
it_lim_it = 0;                          % Iteração máx. estágio de Int.
Tabu_Search_it = 0;                             % Iteração na qual foi pra Diversificação
Intensification_it = 0;
Heuristic_it = 1;


% Duração do banimento
tabu_tenure_best = 0;                   % banimento do melhor movimento