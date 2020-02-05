% Gráfico do Exemplo - Coordenação 2 Relés
clearvars
close all
clc
tic

% CTI, RTC
N_Reles = 2;
CTI = 0.4;
RTC = [300/5, 300/5];

% Faltas máx-mín
FC_max = [3115.0 2010.7];
FC_min = [1741.3 1309.9];

Icc_K1 = FC_max;    % modo máx
Icc_K2 = FC_min;    % modo mín

% Corrente de Carga
n_pontos = 3001;
I_grafico = 200:1:3200;

% Ordem Relés Retaguarda
q = [1; 1];

B1_ret = zeros(n_pontos,1);
B2_pri = zeros(n_pontos,1);

% Curva IEC Inversa, Dial
Ai = [0.14 0.14];
N = [0.02 0.02];
x = [0.1263 0.1000];
% x = [0.2524 0.1000];
% PS = [10 4];
PS = [10 4];

for i = 1:n_pontos
    % Matriz faltas R1
    %B1_pri = Ai(1)/((FC_max(1)/(PS(1)*RTC(1)))^N(1)-1);
    B1_ret(i,1) = Ai(1)/((I_grafico(i)/(PS(1)*RTC(1))).^N(1)-1);
    % Matriz faltas R2
    B2_pri(i,1) = Ai(2)/((I_grafico(i)/(PS(2)*RTC(2))).^N(2)-1);
    %B2_ret = Ai(q(2))/((FC_min(2)/(PS(q(2))*RTC(q(2))))^N(q(2))-1);
end

% Tempos de Atuação

T_pri2 = x(1,2) * B2_pri;
T_ret1 = x(1,1) * B1_ret;
%% Plot

set(groot,'defaultLineLineWidth',1.8)
semilogy(I_grafico, T_pri2)
xlim([0 3115])
ylim([0.1 100])
xlabel('Icc (A)')
ylabel('t (s)')
hold on
grid
semilogy(I_grafico,T_ret1)
hold on
line([2010.7 2010.7],ylim,'color','black','LineWidth',0.8)
line([1309.9 1309.9 ],ylim,'color','black','LineWidth',0.8)
title('Coordenação de 2 Relés de Sobrecorrente')
legend ('Tpri2','Tret1');

t = toc;