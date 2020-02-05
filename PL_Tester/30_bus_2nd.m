%% Sistema de 30 barras - Resolu��o PL

% 68 rel�s - 122 restri��es de seletividade

close all;
clearvars
clc
tic

%% Dados dos rel�s

% Plug Settings 
PS = [3.0 3.5 3.0 2.0 2.5 1.5 3.0 3.5 3.0 4.5 2.0 3.0 2.0 2.5 2.0 2.5 2.5 5.0 3.50 1.50 4.50 5.0 2.5 2.0 5.0 3.5 5.0 3.5 4.5 2.0 2.5 4.0 3.0 1.5 5.0 2.5 4.0 2.0 1.5 5.0 4.5 1.5 4.5 4.0 3.0 3.0 3.5 2.0 3.5 5.0 3.5 2.5 3.5 1.5 2.5 2.5 4.0 4.5 1.5 1.5 3.5 2.0 4.0 3.0 3.0 2.0 4.0 4.5]';

f = zeros(1,136);
CTI = 0.3;
RTC = ones(68,1) * 80;

% Ordem dos rel�s (prim-p / ret-q)
p = [1	2	2	2	3	4	5	5	5	6	6	7	7	7	8	9	9	9	10	10	10	10	11	12	12	13	13	14	14	14	14	15	16	17	17	17	17	18	19	19	19	19	20	21	21	21	21	22	23	24	25	25	25	26	27	27	27	28	29	29	29	30	31	31	31	32	32	33	33	34	35	35	36	36	36	37	37	38	39	40	40	40	41	41	41	42	43	43	43	44	45	46	47	48	49	50	51	52	52	53	53	54	54	55	56	56	57	57	58	59	59	61	62	62	63	63	64	65	65	66	67	68];

q = [4	6	8	10	2	12	1	8	10	11	14	1	6	10	16	1	6	8	13	18	20	22	3	5	14	5	11	9	18	20	22	7	17	9	13	20	22	15	9	13	18	22	24	9	13	18	20	23	19	21	28	30	32	45	26	30	32	49	26	28	32	52	26	28	30	51	54	36	38	40	34	38	39	42	44	34	36	46	33	35	42	44	35	39	44	48	35	39	42	56	37	25	41	50	47	27	29	31	54	31	51	55	58	43	53	58	53	55	62	57	62	57	64	66	61	66	68	61	64	67	63	65];

FC_pr = [9191.74	7600.07	7605.54	4428.08	6239.73	4844.6	6037.99	3537.27	6231.31	4937.31	4607.1	5218.25	5295.19	5414.87	3520.82	3507.76	5118.43	3680.97	5283.83	4304.34	5128.21	3646.18	3709.05	3269.01	9483.26	7821.79	8923.48	5913.93	9469.47	8163	9119.61	7842.95	8445.81	5493.7	9105.57	7593.61	8890.18	6684.02	5426.91	7029.76	7117.7	5467.58	7154.11	5532.86	6593.53	7410.28	5381.07	5430.47	5489.99	5728.76	8101.28	8042.1	7658.65	6330.11	5388.93	6028.42	5865.08	4479	4012.26	13076.1	4445.16	4700.07	4273.66	2642.77	4088.63	2424.48	2495.44	2358.22];

FC_bc = [2645.6	1569.84	1766.69	603.1	1228.33	427.28	1396.38	461.74	1144.4	461.86	1214.07	2454.21	1618.03	1825.06	1055.38	1088.98	2153.96	652.5	2282.09	619.85	2207.8	104.97	691.7	426.35	5051.22	1390.43	3665.21	882.49	4422.86	420.87	2213.48	323.41	2753.35	245.77	3279.28	1070	3355.94	1352.77	870.74	2239.21	2917.77	1203.62	2951.52	1221.22	2032.54	2844.32	2195.58	2234.39	1908.6	3004.33	3384.19	3032.08	2855.42	1629.63	1519.44	2144.43	2040.51	1203.31	9999999	9999999	1430.92	1923.93	1651.43	69.87	1243.7	47.03	854.4	699.51];

Icc_K1 = [FC_pr; FC_bc];

B1 = size(Icc_K1);

% Curva IEC Inversa
Ai = 0.14;
N = 0.02;

options = optimoptions(@linprog,'Display','none');

%% Matrizes Kij

for j = 1:68
    B1(1,j) = Ai/((FC_pr(j)/(PS(j)*RTC(j)))^N-1);  % Rel� prim�rio
    B1(2,j) = Ai/((FC_bc(j)/(PS(j)*RTC(j)))^N-1);  % Rel� retaguarda
end

Kij_geral = zeros(244,136);

for k = 1:244
    if(k <= 122)
        Kij_geral(k,p(k)) = B1(1,p(k));                 % Rprim
        Kij_geral(k,q(k)) = -B1(2,q(k));                % Rret
        else Kij_geral(k,q(k-122)+68) = -1;             % Rtz2
             Kij_geral(k,p(k-122)) = B1(1,p(k-122));    % Rprim 
    end
end

%% Sistema no formato do Linprog

% F.O - Tempo m�nimo de atua��o dos Rel�s (somente Prim�rios)

 for i = 1:136
     if(i <= 68)
         f(1,i) = B1(1,i);
     else f(1,i) = 1;
     end
 end
 
% Valores M�n/M�x de Dial dos Rel�s
lb = ones(1,136);
ub = ones(1,136);

for i = 1:136
    if(i <= 68)
        lb(i) = 0.1;
        ub(i) = 1.1;
    else lb(i) = 0.30;
         ub(i) = 1.50;
    end
end

% Lado direito da Desigualdade
b = ones(244,1);
b = b * -CTI;

A = Kij_geral;

% Exce��es

%% Resolu��o via PL
abc= 1;

[x,fval,exitflag] = linprog(f,A,b,[],[],lb,ub,[],options);

%% Tempos de Atua��o dos Rel�s

T_pri  = zeros(68,1);    % Pr�-aloca��o
T_ret  = zeros(68,1);
T_ret_z2 = zeros(68,1);

for i = 1:68
    T_pri(i,1) = x(i,1) * B1(1,i);
    T_ret(i,1) = x(i,1) * B1(2,i);
    T_ret_z2(i,1) = x(i+68,1);
end

T_coord = ones(122,5);     % Configura��o �nica!

N_violation_OC = 0;
N_violation_tZ2 = 0;

for i = 1:122
    T_coord(i,1) = p(i);
    T_coord(i,2) = q(i);
    T_coord(i,3) = T_pri(p(i));
    T_coord(i,4) = T_ret(q(i));
    T_coord(i,5) = T_ret(q(i)) - T_pri(p(i));
    T_coord(i,6) = T_ret_z2(q(i));
    T_coord(i,7) = T_ret_z2(q(i)) - T_pri(p(i));
    
    if(T_coord(i,5) < 0.30)
        N_violation_OC = N_violation_OC + 1;
    end
        
    if(T_coord(i,7) < 0.30) 
        N_violation_tZ2 = N_violation_tZ2 + 1;
    end
    
end

t = toc;