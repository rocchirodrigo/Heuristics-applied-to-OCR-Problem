%% Tabu Search
clearvars
clc

% Para rodar o programa, informe o código do sistema a ser utilizado,
% em parênteses no campo SYSTEM, e o número de simulações em total.

% (1) Sistema I: 5 Reles, 28 restrições, count_MAX = 11 (8 seletiv. 20 max_min)
% (2) Sistema II: 10 Reles, 58 restrições , count_MAX = 11 (18 seletiv. 40 max_min)
% (3) 3 Barras: 6 Reles, 12 Restrições, count_MAX = 8
% (4) 8 Barras: 14 Reles, 20 Restrições, count_MAX = 7
% (5) 15 Barras: 42 Reles, 82 Restrições, count_MAX = 5
% (6) 2 Relés (Apêndice!)

SYSTEM = 2;
total = 1;
N_it_MAX = 20;

switch SYSTEM
    case 1
        N_Reles = 5;   N_Restricoes = 28; count_MAX = 11;
    case 2
        N_Reles = 10;   N_Restricoes = 58; count_MAX = 11;
    case 3
        N_Reles = 6;   N_Restricoes = 12; count_MAX = 8;
    case 4
        N_Reles = 14;   N_Restricoes = 20; count_MAX = 7;
    case 5
        N_Reles = 42;   N_Restricoes = 82; count_MAX = 5;
    case 6
        N_Reles = 2;   N_Restricoes = 6; count_MAX = 11;
end

exit = zeros(total,1);
PL_total = zeros(total,1);
t_exec = zeros(total,1);
fo_final = zeros(total,1);
PS_inicial = zeros(N_Reles,total);
PS_final = zeros(N_Reles,total);
Dial_final = zeros(N_Reles,total);
evolucao = zeros(N_it_MAX,total);

for vezes = 1:total
    % Limpa todas variáveis exceto as estatísticas.
    clearvars -except SYSTEM t_exec PL_total fo_final PS_final PS_inicial Dial_final total vezes evolucao N_Reles count_MAX N_Restricoes N_it_MAX exit
    clc
    tic
    
    %%  Descomente o Sistema que deseja utilizar
    
    switch SYSTEM
        case 1
            run tabu_radial_i_teste.m;
        case 2
            run tabu_radial_ii_teste.m;
        case 3
            run tabu_3_bus_teste.m;
        case 4
            run tabu_novo_8_bus.m;
%             run tabu_8_bus_teste.m;
        case 5
%             run tabu_novo.m;
            run tabu_15_bus_teste.m;
        case 6
            run tabu_2_reles_teste.m;
    end
    
    % Salva as variáveis referentes àquela execução
    PL_total(vezes,1) = PL_count;
    t_exec(vezes,1) = t;
    fo_final(vezes,1) = best_fo;
%     fo_final(vezes,1) = best_ff;
    
    for i = 1:N_it_MAX
        evolucao(i,vezes) = it_results(i,1);
    end
    
    for i = 1:N_Reles
        PS_final(i,vezes) = best_move(i,1);
        Dial_final(i,vezes) = best_x(i,1);
    end
    
end

%% Estatísticas para TCC (100 Execuções):

% Função Objetivo
FO_min = min(fo_final);
FO_media = mean(fo_final);
FO_max = max(fo_final);
FO_desvio = std(fo_final);

% PL Solucionados
PL_min = min(PL_total);
PL_media = mean(PL_total);
PL_max = max(PL_total);
PL_desvio = std(PL_total);

% Tempo de Execução
T_Exe_min = min(t_exec);
T_Exe_media = mean(t_exec);
T_Exe_max = max(t_exec);
T_Exe_desvio = std(t_exec);