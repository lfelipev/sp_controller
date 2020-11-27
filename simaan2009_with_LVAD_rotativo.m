clearvars -except SP_constant w_rpm_constant EDP_constant COvec_constant SP_controlador w_rpm_controlador EDP_controlador COvec_controlador
% close all
clc

% Simulation Time;
start_t = 0;
passo   = 0.0001;
end_t   = 60;

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
Rm  = 0.1; % Rm-mitral valve open;(mmHg.sec/mL)
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
alfa = -3.5;

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
EDP(1) = Pve(1);
Aux(1) = 40;
estado(1) = 0;
enable_preload = 1;

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
        if i < 80000 % 1a constante
            Rm = 0.1;
            cae(i+1) = Cae;
            rm(i+1) = Rm;
        % 8000 -> 12000 // Rm = Variavel e Cae = variavel
        elseif i >= 80000 && i < 120000 % rampa de subida
            Rm = -2.375e-6*i+0.29000000000000004;
            %Cae = 1.8e-3*i -136;
            cae(i+1) = Cae;
            rm(i+1) = Rm;
        % 12000 -> 32000 // Rm = 0.001 e Cae = 600
        elseif i >= 120000 && i < 320000 % 2a constante
            Rm = 0.005;
            %Cae = 80;
            cae(i+1) = Cae;
            rm(i+1) = Rm;
        % 32000 -> 40000 // Rm = variavel e Cae = 600
        elseif i >= 320000 && i <= 400000 % rampa de descida
            Rm = 3.0625e-6*i - 0.9749999999999999;
            %Cae = 80;
            cae(i+1) = Cae;
            rm(i+1) = Rm;
        elseif i >= 400000
            Rm = 0.25;
            %Cae = 80;
            cae(i+1) = Cae;
            rm(i+1) = Rm;
        end
    end
   
    % Enchimento - 1
    % Dm = 1, Da = 0
    
    % Aqui é calculado o End Diastolic Volume
    % Quando o estado anterior é enchimento e o estado atual é CI
    
    % Contração Isovolumétrica - 2
    % Dm = 0, Da = 0
    
    % Ejeção - 3
    % Dm =0, Da = 1
    
    % Relaxamento Isovolumétrico - 4
    % Dm = 0, Da = 0
    
    w = (w_rpm(i)*2*pi/60); % rad/s
    % w = (w_rpm/60)*2*pi; % rad/s

    %Rk
    if Pve(i) > 1 % x1_
        Rk = 0;
    else
        Rk = 0;
        %Rk = alfa*(Pve(i) - 1);
    end
    RR = Ri + Ro + Rk + Bo;
    %RR = Ri + Ro + Bo;
    
    % Matrix A
    % lambda = E(i) - Emin*En(i);
    a13 = (Da/Ra)*(E(i));
    a33 = -((Dm/Rm)+(Da/Ra))*E(i);
    a53 = (Dm/Rm)*E(i);
    a55 = -(1/Rs + Dm/Rm);

    %    x = [   x1     x2    x3    x4     x5     x6  ]
    %    x = [   Pao    Qa    Vve   Pas    Pae   Qvad ]
    A(1,:) = [ -Da/Ra   -1    a13    0      0      1  ]/Cao;
    A(2,:) = [    1     -Rc    0    -1      0      0  ]/Ls;
    A(3,:) = [  Da/Ra    0    a33    0    Dm/Rm   -1  ];
    A(4,:) = [    0      1     0   -1/Rs   1/Rs    0  ]/Cs;
    A(5,:) = [    0      0    a53   1/Rs   a55     0  ]/Cae;
    A(6,:) = [   -1      0    E(i)   0      0    -RR  ]/LL;

    % Matrix B
    B(1,1) = (-(Da/Ra)*E(i)*Vo)/Cao;
    B(2,1) = 0;
    B(3,1) = (Dm/Rm + Da/Ra)*E(i)*Vo;
    B(4,1) = 0;
    B(5,1) = (-(Dm/Rm)*E(i)*Vo)/Cae;
    B(6,1) = (-(E(i)*Vo) - B2*w^2)/LL;

    % O fator w^2 está correto?
    x6dot(i) = (-1*x(1) + E(i)*x(3) -RR*x(6) + (-(E(i)*Vo) - B2*w^2))/LL;

    PIP(i) = Pve(i) - Li*x6dot(i) - Ri*x(6);
    if i > 4
        PIPft(i) = 2.9987*PIPft(i-1) - 2.9987*PIPft(i-2) + 0.9987*PIPft(i-3) + 0.06e-8*PIP(i);
    else
        PIPft(i) = PIP(i);
    end
    PIPf(i) = PIP(i);
    
    if i > 1
        d_PIPf(i) = (PIPf(i) - PIPf(i-1))/passo;
        if d_PIPf(i-1) >= 0.01 && d_PIPf(i) < -0.01
            SP(i) = PIPf(i);
        else
            SP(i) = SP(i-1);
        end
    end
    
    if i > 1
        if(estado_anterior ~= estado_atual && estado_anterior == 2)
            % Encontrar a EDP
            EDP(i) = PIP(i);
            Aux(i) = 40;
        else
            EDP(i) = EDP(i-1);
            Aux(i) = 30;
        end
    end

    x = runkut4(passo,A,x,B);
    
    %w_rpm(i+1) = 12000;
    %w_rpm(i+1) = 12000 + 100*T(i);
    %w_rpm(i+1) = 100*passo + w_rpm(i); % rpm
    %w_rpm(i+1) = 8195;
    phase = 0; % pi/4;
    
    SPref = 87.1316;
    
    w_rpm(i+1) = 8000;
    
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
Aux(i+1) = Aux(i);

% PIP plot
figure(1)
plot(T, SP, '-k','LineWidth', 2);
hold on
plot(T,PIP, 'Color', '#777777')
ylim([0 150])
xlim([0 45])
xticks([0 20 40 60])
yticks([0 50 100 150])
legend('Pump inlet pressure', 'Systolic pressure')
xlabel('Time (s)')
ylabel('PIP (mm Hg)')
title('Preload variation')

% Cardiac Output Plot
figure(2)
subplot(2, 1, 1)
plot(T, COvec, 'Color', '#777777','LineWidth', 2)
ylim([0 10])
yticks([0 2.5 5 7.5 10])
ylabel('CO (L/min)')
xticks([0 20 40 60])
grid on
xlim([0 45])
title('Preload variation')
legend('k_sp = 40 rpm/mm Hg', 'Orientation','horizontal')
legend('boxoff')

subplot(2,1,2)
plot(T, SP,'Color','#777777','LineWidth', 2)
ylim([0 150])
yticks([0 50 100 150])
xticks([0 20 40 60])
grid on
xlim([0 45])
xlabel('Time(s)')
ylabel('SP (mm Hg)')
legend('k_sp = 40 rpm/mm Hg', 'Orientation','horizontal')
legend('boxoff')

%%
[PIP_phis, SP_phis, COvec_phis, EDP_phis] = physiological_simaan(1, 60);
[PIP_controlador, SP_controlador, COvec_controlador, w_rpm_controlador, EDP_controlador] = sp_controller_simaan(1, 60);
COvec_constant = COvec;
SP_constant = SP;
w_rpm_constant = w_rpm;
EDP_constant = EDP;

%% Plot do petrou
figure(2)

subplot(4,1,1);
plot(T, w_rpm_controlador/1000, 'k')
hold on
plot(T, w_rpm_constant/1000, '-k', 'LineWidth', 2)
title('Preload variation')
ylabel('Pump Speed (krpm)')
ylim([6 9.5])
xlim([0 45])
yticks([2 3 4 5 6])
xticks([0 20 40 60])
legend('SP controller', 'Constant speed');
grid on

subplot(4,1,2);
plot(T,COvec_controlador, 'k')
hold on
plot(T,COvec_constant, '-k', 'LineWidth', 2);
ylabel('CO (L/min)')
hold on
plot(T,COvec_phis, ':k', 'LineWidth', 2);
ylim([1 9])
xlim([0 45])
yticks([1 3 5 7 9])
xticks([0 20 40 60])
legend('SP controller', 'Constant speed', 'Physiological');
grid on

subplot(4,1,3);
plot(T,SP_controlador, 'k')
hold on
plot(T, SP_constant, '-k','LineWidth', 2);
hold on
plot(T,SP_phis, ':k', 'LineWidth', 2);
ylabel('SP (mm Hg)')
ylim([50 150])
xlim([0 45])
yticks([0 50 100 150 200])
xticks([0 20 40 60])
legend('SP controller', 'Constant speed', 'Physiological');
grid on

subplot(4,1,4);
plot(T,EDP_controlador, 'k')
hold on
plot(T, EDP_constant, '-k', 'LineWidth', 2);
hold on
plot(T,EDP_phis, ':k', 'LineWidth', 2);
ylabel('EDP (mm Hg)')
xlabel('Time (s)')
legend('SP controller', 'Constant speed', 'Physiological');
ylim([30 90])
xlim([0 45])
xticks([0 20 40 60])
yticks([0 10 20 30 40 50 60 70 80 90])
grid on
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