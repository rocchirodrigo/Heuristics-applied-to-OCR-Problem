% Gráfico n.Iterações x FO
clearvars
close all
clc

iter = 1:1:10;
FO = [13.3913 12.2892 12.2244 12.2149 12.2149 12.2149 12.2149 12.2149 12.2149 12.2149];
%% Plot

set(groot,'defaultLineLineWidth',1.5)
plot(iter,FO)
xlim([1 10])
 ylim([12 14])
xlabel('Iterações')
ylabel('FO (s)')
hold on
grid
legend ('FO');