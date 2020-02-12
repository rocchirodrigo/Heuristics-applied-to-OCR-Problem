%% Sistema de 3 barras - Resolução PL

close all;
clearvars
clc
tic

%% Dados dos relés

N_Reles = 5;
N_Restricoes = 28;
count_MAX = 11;

% Plug Settings
PS = [5.50 4 5.50 4 4];
PS_disp = zeros(count_MAX,1);

for i = 1:count_MAX
    PS_disp(i,1) = 2.5 + 0.75*(i-1);
end


CTI = 0.4;
RTC = [300/5, 300/5, 100/5, 200/5, 100/5];

% Ordem dos relês (prim-p / ret-q)
p = [1; 2; 3; 4; 5];
q = [1; 1; 1; 2; 3];

% Odem dos relês primários na matriz
r = [1; 2; 3; 4];
s = [2; 3; 4; 5];

% Odem dos relês retaguarda na matriz
t = [1; 2; 3; 4];
u = [1; 1; 2; 3];

% Faltas máx-mín
FC_max = [3115.0 2010.7 2010.7 1512.5 878.4];
FC_min = [1741.3 1309.9 760.7 500.3 325.1];

Icc_K1 = FC_max;    % modo máx
Icc_K2 = FC_min;    % modo mín

% Corrente de Carga
I_load = [199.5 130.8 68.7 100.7 50.0];

% Limites PS - Carga
FC = 1.50;
lim_inf = ones(N_Reles,1);
lim_sup = ones(N_Reles,1);

% Restrição do valor de PS
for i = 1:N_Reles
    lim_inf(i,1) = (FC*I_load(i))/RTC(i);
    lim_sup(i,1) = FC_min(i)/RTC(i);
end

% Valores mínimos permitidos pelo Relé
PS_min = interp1(PS_disp, PS_disp, lim_inf,'next');
PS_max = interp1(PS_disp, PS_disp, lim_sup,'nearest','extrap');
rest_carga = ones(11,N_Reles);
rest_carga2 = ones(11,N_Reles);

for i = 1:11
    for j = 1:N_Reles
        rest_carga(i,j) = PS_disp(i,1);
        rest_carga2(i,j) = PS_disp(i,1);
        if((rest_carga(i,j)) < lim_inf(j,1) || (rest_carga(i,j)) > lim_sup(j,1))
            rest_carga(i,j) = 1;    % Fora da restrição
        else rest_carga(i,j) = 0;   % Cumpre a restrição
        end
    end
end

B1 = size(Icc_K1);
B2 = size(Icc_K2);
f = zeros(1,N_Reles);

% Curva IEC Inversa
Ai = [13.50 0.14 0.14 13.50 80.00];
Ni = [1.00 0.02 0.02 1.00 2.00];

options = optimoptions(@linprog,'Display','none');

%% Matrizes Kij

% Matrizes Kij
for j = 1:N_Reles
    % Matriz faltas max
    B1(1,j) = Ai(j)/((FC_max(j)/(PS(j)*RTC(j)))^Ni(j)-1);
    B1(2,j) = Ai(q(j))/((FC_max(j)/(PS(q(j))*RTC(q(j))))^Ni(q(j))-1);
    % Matriz faltas min
    B2(1,j) = Ai(j)/((FC_min(j)/(PS(j)*RTC(j)))^Ni(j)-1);
    B2(2,j) = Ai(q(j))/((FC_min(j)/(PS(q(j))*RTC(q(j))))^Ni(q(j))-1);
end

Kij_geral = zeros(N_Restricoes,N_Reles);

for k = 1:4
    Kij_geral(k,s(k)) = B1(1,s(k)); % prim. max
    Kij_geral(k,u(k)) = -B1(2,s(k)); % ret. max
    Kij_geral(4+k,s(k)) = B2(1,s(k)); % prim. min
    Kij_geral(4+k,u(k)) = -B2(2,s(k)); % ret. min
end

for k = 1:5
    % Tmin
    Kij_geral(8+k,k) = -B1(1,k);
    Kij_geral(13+k,k) = -B2(1,k);
    % Tmax
    Kij_geral(18+k,k) = B1(1,k);
    Kij_geral(23+k,k) = B2(1,k);
end

%% Sistema no formato do Linprog

% F.O - Tempo mínimo de atuação dos Relés (somente Primários)

for i = 1:N_Reles
    f(1,i) = B1(1,i) + B2(1,i);
end

% Valores Mín/Máx de Dial dos Relés R1-R5
lb = ones(1,N_Reles) * 0.1;
ub = ones(1,N_Reles) * 10;

% Lado direito da Desigualdade
b = ones(N_Restricoes,1);

% CTI
for i = 1:8
    b(i,1) = b(i,1) * -0.4;
end

%T min
for i = 1:10
    b(i+8,1) = b(i+8,1) * -0.05;
end

% Tmax
for i = 1:10
    b(i+18,1) = b(i+18,1) * 2;
end

A = Kij_geral;

%% Resolução via PL

[x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);

%% Tempos de Atuação dos Relês

T_pri1 = zeros(N_Reles,1);    % Pré-alocação
T_ret1 = zeros(N_Reles,1);
T_pri2 = zeros(N_Reles,1);
T_ret2 = zeros(N_Reles,1);

for i = 1:N_Reles;
    T_pri1(i,1) = x(i,1) * B1(1,i);
    T_ret1(i,1) = x(q(i),1) * B1(2,i);
    T_pri2(i,1) = x(i,1) * B2(1,i);
    T_ret2(i,1) = x(q(i),1) * B2(2,i);
end

T_coord1 = zeros(N_Reles,1);     % Configuração única!
T_coord2 = zeros(N_Reles,1);
violation1 = zeros(N_Restricoes,3);
violation2 = zeros(N_Restricoes,3);
index = 1;

for i = 1:N_Reles
    T_coord1(i,1) = p(i);
    T_coord1(i,2) = q(i);
    T_coord1(i,3) = T_pri1(i);
    T_coord1(i,4) = T_ret1(i);
    T_coord1(i,5) = T_coord1(i,4) - T_coord1(i,3);
    
    if(T_coord1(i,5) < CTI)  % Violações na Coordenação
        violation1(index,1) = p(i);
        violation1(index,2) = q(i);
        violation1(index,3) = T_coord1(i,5);
        index = index+1;
    end
end

for i = 1:N_Reles
    T_coord2(i,1) = p(i);
    T_coord2(i,2) = q(i);
    T_coord2(i,3) = T_pri2(i);
    T_coord2(i,4) = T_ret2(i);
    T_coord2(i,5) = T_coord2(i,4) - T_coord2(i,3);
    
    if(T_coord2(i,5) < CTI)  % Violações na Coordenação
        violation2(index,1) = p(i);
        violation2(index,2) = q(i);
        violation2(index,3) = T_coord2(i,5);
        index = index+1;
    end
end

t = toc;