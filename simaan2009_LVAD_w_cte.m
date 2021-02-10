clear
close all
clc

global Da Dm Rs Rm Ra Rc Cae Cs Cao Ls RR LL Vo B2 alpha

%% Simulation Time;
start_t = 0;
passo   = 0.0001;
end_t   = 120;

% Uses the already created Time scale
T = start_t:passo:end_t;
n = length(T);

%% Cardiovascular system
HR = 90;
Emax = 1.50;
Emin = 0.06;
En = Elastance(T,passo,HR,end_t);
E = (Emax - Emin)*En + Emin;

% Cardiovascular system model parameters (from Simaan2009);
Rs  = 1.0000; % (0.83-normal,weak; 1.4-severly weak without pump; 0.83-severly weak with pump)(mmHg.sec/mL)
Rm  = 0.1000; % Rm-mitral valve open;(mmHg.sec/mL)
Ra  = 0.0010; % Ra-aortic valve open;(mmHg.sec/mL)
Rc  = 0.0398; % Rc-characteristic resistance;(mmHg.sec/mL)
Cae = 4.4000; % Cr-pulmonary compliance;(mL/mmHg)
Cs  = 1.3300; % Systemic Complinace (ml/mmHg)
Cao = 0.0800; % Aortic Complinace (ml/mmHg)
Ls  = 0.0005; % Ls-inertance of blood in aorta;(mmHg.sec^2/mL)

% LVAD parameters
Ri = 0.0677;
Ro = 0.0677;
Li = 0.0127;
Lo = 0.0127;

% Pump
Bo = 0.17070;
B1 = 0.02177;
B2 = -9.9025e-5;

LL = Li + Lo + B1;
alpha = -3.5;

Vo = 10;

% Initial Conditions
Pao(1)  = 90;
Qa(1)   = 0;
Vve(1)  = 140;
Pas(1)  = 90;
Pae(1)  = 5;
Qvad(1) = 0; % x6 - LVAD flow

%x = [  x1     x2      x3      x4      x5        x6  ]';
x =  [Pao(1)  Qa(1)  Vve(1)  Pas(1)  Pae(1)   Qvad(1)]';

% Initial states of diodes
Dm = 0; Da = 0;

CO = 0;
EDV = 0;
ESV = 0;

% constant speed

w_rpm = 7515;
w = (w_rpm*2*pi/60);

estado_atual = 3;
SP      = zeros(1,length(T));

rs  = ones(1, length(T));
rm  = zeros(1, length(T));
cae = zeros(1, length(T));

Pao = zeros(1, length(T));
Qa  = zeros(1, length(T));
Vve = zeros(1, length(T));
Pas = zeros(1, length(T));
Pae = zeros(1, length(T));
Qvad = zeros(1, length(T));

Pve  = zeros(1, length(T));

PIP = zeros(1, length(T));
d_PIP  = zeros(1, length(T));
COvec  = zeros(1, length(T));
EDVvec = zeros(1, length(T));
ESVvec = zeros(1, length(T));
estado = zeros(1, length(T));

estado(1) = 0;
enable_preload = 1;
rm(1) = 0.1;
SP(1) = 1;
PIP(1) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Pae(i) >= Pve(i)
    Dm = 1;
else
    Dm = 0;
end

if Pve(i) >= Pao(i)
    Da = 1;
else
    Da = 0;
end

estado_anterior = estado_atual;

if (Dm == 1 && Da == 0) && (estado_atual == 4)
    estado_atual = 1; % se o estado atual é Relaxamento(4), vai p/ Enchimento(1)
    ESV = Vve(i);
end
if (Dm == 0 && Da ==0) && (estado_atual == 1)
    estado_atual = 2; % se o estado atual é Enchimento(1), vai p/ Contração Isovolumétrica(2)
end
if (Dm == 0 && Da == 1) && (estado_atual == 2)
    estado_atual = 3; % Se o estado atual é Contração Isovolumétrica(2), vai p/ Ejeção(3)
    EDV = Vve(i);
end
if (Dm == 0 && Da == 0) && (estado_atual == 3)
    estado_atual = 4; % Se o estado atual é Ejeção(3), vai p/ Relaxamento Isovolumétrico(4)
end

COvec(i+1) = ((EDV-ESV) * HR) / 1000;
EDVvec(i+1) = EDV;
ESVvec(i+1) = ESV;

estado(i+1) = estado_atual;

% Variação da Pré-carga
if enable_preload
    if i < 160000 % 1a constante
        Rm = 0.1;
        cae(i+1) = Cae;
        rm(i+1) = Rm;
    % 8000 -> 12000 // Rm = Variavel e Cae = variavel
    elseif i >= 160000 && i < 200000 % rampa de subida
        Rm = -5e-7*i+0.18;
        %Cae = 1.8e-3*i -136;
        cae(i+1) = Cae;
        rm(i+1) = Rm;
    elseif i>= 200000 && i < 300000
        Rm = 0.08;
        rm(i+1) = Rm;
    elseif i >= 300000 && i < 400000
        Rm = 1e-7*i+0.05;
        rm(i+1) = Rm;
    elseif i>= 400000
        Rm = 0.09;
        rm(i+1) = Rm;
    end
end

% Variação da Pós-carga
if i >= 600000 && i < 700000
    rs(i) = -5e-6 * i + 4;
elseif i>= 700000 && i < 800000
    rs(i) = 0.5;
elseif i >= 800000 && i<900000
    rs(i) = 2.5e-6 * i - 1.5;
elseif i >= 900000
    rs(i) = 0.75;
end
Rs = rs(i);

% Rk
if Pve(i) > 1 % x1_
    Rk = 0;
else
    Rk = 0; % alpha*(Pve(i) - 1);
end
RR = Ri + Ro + Rk + Bo;

xdot = xdot_fun(x,E(i),w);

PIP(i) = Pve(i) - Li*xdot(6,1) - Ri*x(6);
if i > 1
    d_PIP(i) = (PIP(i) - PIP(i-1))/passo;
    if d_PIP(i-1) >= 0.01 && d_PIP(i) < -0.01
        SP(i) = PIP(i);
    else
        SP(i) = SP(i-1);
    end
end

x = runkut42(x,xdot,E(i),w,passo);

Pao(i+1) = x(1);
Qa(i+1)  = x(2);
Vve(i+1) = x(3);
Pas(i+1) = x(4);
Pae(i+1) = x(5);
Qvad(i+1)= x(6);

Pve(i+1) = E(i+1)*(Vve(i+1) - Vo);
end
PIP(i+1) = PIP(i);
rs(i+1) = 0.75;
SP(i+1) = SP(i);

COvec_w_cte = COvec;
SP_w_cte    = SP;

%%
save('simaan2009_LVAD_w_cte.mat','COvec_w_cte','SP_w_cte')

%% PIP plot
figure(1)
plot(T, SP, '-k','LineWidth', 2);
hold on
plot(T,PIP, 'color',[0.5 0.5 0.5]);
ylim([0 120])
xlim([0 120])
xticks([0 20 30 40 60 70 80 90 120])
yticks([0 50 100 120])
legend('PIP', 'SP')
xlabel('Time (s)','interpreter','latex')
ylabel('Pressure (mmHg)','interpreter','latex')
set(gca,'FontSize',14)
set(gca,'fontname','times')
grid on

%%
figure(2)
subplot(2, 1, 1);
plot(T, rm,'-k','LineWidth', 2)
grid on
ylim([0.07 0.11])
yticks([0.08 0.09 0.1])
xticks([0 16 20 30 40 120])
ylabel('$R_m$','interpreter','latex')
set(gca,'FontSize',16)
set(gca,'fontname','times')

subplot(2, 1, 2);
plot(T, rs, '-k','LineWidth', 2)
grid on
ylim([0.4 1.1])
xticks([0 60 70 80 90 120])
yticks([0.5 0.75 1])
xlabel('Time (s)')
ylabel('$R_s$','interpreter','latex')
set(gca,'FontSize',16)
set(gca,'fontname','times')
%%
% Cardiac Output Plot
figure(3)
subplot(2, 1, 1)
plot(T, COvec, 'Color', [0.5 0.5 0.5],'LineWidth', 2)
ylabel('CO (L/min)')
xticks([0 20 40 60 100])
grid on
axis([0 120 3 6])
title('Preload variation')
legend('k_sp = 40 rpm/mm Hg', 'Orientation','horizontal')
legend('boxoff')

subplot(2,1,2)
plot(T, SP,'Color',[0.5 0.5 0.5],'LineWidth', 2)
ylim([0 150])
yticks([0 50 100 150])
xticks([0 20 40 60])
grid on
xlim([0 100])
xlabel('Time(s)')
ylabel('SP (mm Hg)')
legend('k_sp = 40 rpm/mm Hg', 'Orientation','horizontal')
legend('boxoff')