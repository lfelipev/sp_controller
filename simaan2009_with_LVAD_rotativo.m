clearvars -except w_rpm_2 w_rpm_constant EDP_constant COvec_constant SP_controlador...
           w_rpm_controlador EDP_controlador COvec_controlador
% close all
clc

global Da Dm Rs Rm Ra Rc Cae Cs Cao Ls RR LL Vo B2 alpha

% Simulation Time;
start_t = 0;
passo   = 0.0001;
end_t   = 40;

%Uses the already created Time scale
T = start_t:passo:end_t;
n = length(T);

% Cardiovascular system
HR = 75; %75
Emax = 2; %1.5
Emin = 0.05;
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
B1 = 0.02177;6100
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

Pve(1) = E(1)*(Vve(1) - Vo);

%x = [  x1     x2      x3      x4      x5        x6  ]';
x =  [Pao(1)  Qa(1)  Vve(1)  Pas(1)  Pae(1)   Qvad(1)]';

% Initial states of diodes
Dm = 0; Da = 0;

Da_(1) = Da;
Dm_(1) = Dm;
CO = 0;
EDV = 0;
ESV = 0;
% Simulation
w_rpm = zeros(1, end_t/passo+1);
w_rpm(1) = 6100;
estado_atual = 3;
d_PIP(1) = 0; % derivada do PIP
SP      = zeros(1,length(T));
rs = ones(1, length(T));
EDP(1) = Pve(1);
Aux(1) = 40;
estado(1) = 0;
enable_preload = 1;
rm(1) = 0.1;
SP(1) = 1;
PIP(1) = 0;

for i = 1:n-1

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
%SV = EDV - ESV;
estado(i+1) = estado_atual;

if enable_preload
    % Varia?ão da Pré-carga
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

w = (w_rpm(i)*2*pi/60); % rad/s

%Rk
if Pve(i) > 1 % x1_
    Rk = 0;
else
    Rk = 0;
    %Rk = alpha*(Pve(i) - 1);
end
RR = Ri + Ro + Rk + Bo;
%RR = Ri + Ro + Bo;

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

if i > 1
    if(estado_anterior ~= estado_atual && estado_anterior == 2)
        % Encontrar a EDP
        EDP(i) = PIP(i);
    else
        EDP(i) = EDP(i-1);
    end
end

x = runkut42(x,xdot,E(i),w,passo);

phase = 0; % pi/4;

w_rpm(i+1) = 8660;

Pao(i+1) = x(1);
Qa(i+1)  = x(2);
Vve(i+1) = x(3);
Pas(i+1) = x(4);
Pae(i+1) = x(5);
Qvad(i+1)= x(6);

Da_(i+1) = Da;
Dm_(i+1) = Dm;

Pve(i+1) = E(i+1)*(Vve(i+1) - Vo);
end
EDP(i+1) = EDP(i);
PIP(i+1) = PIP(i);
rs(i+1) = 0.75;
SP(i+1) = SP(i);

% PIP plot
% figure(1)
% plot(T, SP, '-k','LineWidth', 2);
% hold on
% % plot(T,PIP, 'Color', '#777777');
% plot(T,PIP, 'yellow');
% ylim([0 120])
% xlim([0 120])
% xticks([0 20 30 40 60 70 80 90 120])
% yticks([0 50 100 120])
% legend('PIP', 'SP')
% xlabel('Time (s)','interpreter','latex')
% ylabel('Pressure (mmHg)','interpreter','latex')
% set(gca,'FontSize',14)
% set(gca,'fontname','times')
% grid on

%%
figure(1)
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
figure(2)
subplot(2, 1, 1)
plot(T, COvec, 'Color', 'yellow','LineWidth', 2)
ylim([0 10])
yticks([0 2.5 5 7.5 10])
ylabel('CO (L/min)')
xticks([0 20 40 60 100])
grid on
xlim([0 100])
title('Preload variation')
legend('k_sp = 40 rpm/mm Hg', 'Orientation','horizontal')
legend('boxoff')
pause
% 
% subplot(2,1,2)
% plot(T, SP,'Color','yellow','LineWidth', 2)
% ylim([0 150])
% yticks([0 50 100 150])
% xticks([0 20 40 60])
% grid on
% xlim([0 100])
% xlabel('Time(s)')
% ylabel('SP (mm Hg)')
% legend('k_sp = 40 rpm/mm Hg', 'Orientation','horizontal')
% legend('boxoff')


%%
% COphy + Petrou + Constant
% tempo = 120;
% [PIP_phis, SP_phis, COvec_phis, EDP_phis, Pas_phis, rs_phis] = physiological_simaan(1, tempo);
% 
% [rs, Qvad, Vve, EDVvec, ESVvec, Pve_controlador, x6dot, PIP_controlador,...
%     SP_controlador, COvec_controlador, w_rpm_controlador, EDP_controlador,...
%     Pas_controlador] = sp_controller_simaan(1, tempo);

%[kspvec, Pve_vgcontrolador, PIP_vgcontrolador, SP_vgcontrolador,...
%     COvec_vgcontrolador, w_rpm_vgcontrolador, EDP_vgcontrolador,...
%     Pas_vgcontrolador] = vgsp_controller_simaan(1, tempo, COvec_phis, SP_phis);

COvec_constant = COvec;
SP_constant = SP;
w_rpm_constant = w_rpm;
EDP_constant = EDP;
Pas_constant = Pas;

%%
% plot(T, COvec_phis, 'Color', '#777777', 'LineWidth', 3)
% hold on
% plot(T, COvec_vgcontrolador, ':k', 'LineWidth', 2)
% hold on
% legend('CO_{phy}', 'CO', 'Location', 'southeast')
% ylim([3 6])
% ylabel('CO (L/min)', 'interpreter','latex')
% xlabel('Time(s)', 'interpreter','latex')
% set(gca,'FontSize',14)
% set(gca,'fontname','times')
% grid on
%%
% figure(3)
% plot(T, COvec_phis, 'Color', 'yellow', 'LineWidth', 3)
% hold on
% plot(T, COvec_controlador, ':k', 'LineWidth', 2)
% hold on
% plot(T, COvec_constant, '--k', 'LineWidth', 2)
% legend('CO_{phy}', 'CO_{SP}', 'CO_{con}', 'Location', 'southeast')
% ylim([3 6])
% ylabel('CO (L/min)', 'interpreter','latex')
% xlabel('Time(s)', 'interpreter','latex')
% set(gca,'FontSize',14)
% set(gca,'fontname','times')
% grid on

%% Plot do petrou
% figure(1)
% 
% subplot(2,1,1);
% plot(T, w_rpm_controlador/1000, 'k')
% hold on
% plot(T, w_rpm_constant/1000, '-k', 'LineWidth', 2)
% title('Preload variation')
% ylabel('Pump Speed (krpm)')
% ylim([6 9.5])
% xlim([0 100])
% yticks([2 3 4 5 6])
% xticks([0 20 40 60])
% legend('SP controller', 'Constant speed');
% grid on
% 
% subplot(2,1,2);
% plot(T,COvec_controlador, 'k')
% hold on
% plot(T,COvec_constant, '-k', 'LineWidth', 2);
% ylabel('CO (L/min)')
% hold on
% plot(T,COvec_phis, ':k', 'LineWidth', 2);
% ylim([0 9])
% xlim([0 100])
% yticks([1 3 5 7 9])
% xticks([0 20 40 60])
% legend('SP controller', 'Constant speed', 'Physiological');
% grid on
% hold on
% plot(T, rs)
% 
% figure(2)
% subplot(2,1,1);
% plot(T,SP_controlador, 'k')
% hold on
% plot(T, SP_constant, '-k','LineWidth', 2);
% hold on
% plot(T,SP_phis, ':k', 'LineWidth', 2);
% ylabel('SP (mm Hg)')
% ylim([50 150])
% xlim([0 100])
% yticks([0 50 100 150 200])
% xticks([0 20 40 60])
% legend('SP controller', 'Constant speed', 'Physiological');
% grid on
% 
% subplot(2,1,2);
% plot(T,EDP_controlador, 'k')
% hold on
% plot(T, EDP_constant, '-k', 'LineWidth', 2);
% hold on
% plot(T,EDP_phis, ':k', 'LineWidth', 2);
% ylabel('EDP (mm Hg)')
% xlabel('Time (s)')
% legend('SP controller', 'Constant speed', 'Physiological');
% ylim([30 90])
% xlim([0 100])
% xticks([0 20 40 60])
% yticks([0 10 20 30 40 50 60 70 80 90])
% grid on
% 
% figure(3)
% subplot(3, 1, 1);
% plot(T, Pas_phis, 'k');
% legend('Phisiological');
% ylabel('Pas (mm Hg)')
% grid on
% title('Pressão Arterial Sistêmica')
% 
% subplot(3, 1, 2);
% plot(T, Pas_controlador, 'k')
% ylabel('Pas (mm Hg)')
% legend('Constant speed');
% grid on
% 
% subplot(3, 1, 3);
% plot(T, Pas_constant, 'k');
% ylabel('Pas (mm Hg)')
% xlabel('Time (s)')
% grid on
% legend('SP controller');


%% Filtering
% 
% fs = 1/passo;
% S = PIP;
% n = length(S);
% X = fft(S);
% f = (0:n-1)*(fs/n);     %frequency range
% power = abs(X).^2/n;    %power
% 
% Y = fftshift(X);
% fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
% powershift = abs(Y).^2/n;     % zero-centered power
% % plot(fshift,powershift)


%%
% plot(PIPf)
% fs = 1/passo;
% PIPf = fft(PIP);
% n = length(PIP);
% f = (0:n-1)*(fs/n);
% power = abs(PIPf).^2/n;
% plot(f, power)
% fc = 1/(2*pi*R*C)
% clc
% Filtro passa baixa para o PIP
%PIPf = lowpass(PIP, 1e-5);
% R = 5e5;
% C = 10;
% for i=1:n
%     PIPf(i) = PIP(i)*(1-exp(-i/(R*C)));
% end

% plot(T, PIPf)


%%
% SP      = zeros(1,length(T));
% threshold = 0;
% teste = 0;

% Calcula a primeira e segunda derivada do PIP filtrado
% for i=1:n
%     if i > 1
%         d_PIP(i) = (PIPf(i)-PIPf(i-1))/passo;
%         dd_PIP(i) = (d_PIP(i)-d_PIP(i-1))/passo;     
%     
%     if dd_PIP(i-1)/1e5 >= -16 && dd_PIP(i)/1e5 <= -16 
%         SP(i) = PIPf(i);
%     else
%        SP(i) = SP(i-1); 
%     end
%    end
% end
% plot(T, dd_PIP/1e5, T, PIPf, T, SP)
% grid on