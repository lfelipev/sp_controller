clear all
% close all
clc

%% Simulation Time;
start_t = 0;
passo   = 0.0001;
end_t   = 20;

%Uses the already created Time scale
T = start_t:passo:end_t;
n = length(T);

%% Cardiovascular system
HR = 60;
Emax = 2.5;
Emin = 0.06;
En = Elastance(T,passo,HR,end_t);
E = (Emax - Emin)*En + Emin;

% Cardiovascular system model parameters (from Simaan2009);
Rs  = 1.0000; % (0.83-normal,weak; 1.4-severly weak without pump; 0.83-severly weak with pump)(mmHg.sec/mL)
%Rm  = 0.0050; % Rm-mitral valve open;(mmHg.sec/mL)
Ra  = 0.0010; % Ra-aortic valve open;(mmHg.sec/mL)
Rc  = 0.0398; % Rc-characteristic resistance;(mmHg.sec/mL)
% Cae = 4.4000; % Cr-pulmonary compliance;(mL/mmHg)
Cs  = 1.3300; % Systemic Complinace (ml/mmHg)
Cao = 0.0800; % Aortic Complinace (ml/mmHg)
Ls  = 0.0005; % Ls-inertance of blood in aorta;(mmHg.sec^2/mL)

%Cae = 3.50;

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

Vo = 20;

% Initial Conditions
Pao(1)  = 90;
Qa(1)   = 0;
Vve(1)  = 140;
Pas(1)  = 90;
Pae(1)  = 5;
Qvad(1) = 0; % x6 - LVAD flow

Pve(1) = E(1)*(Vve(1) - Vo);

w_rpm = zeros(1,end_t/passo);
w_rpm(1) = 12000;
% w_rpm = 12000;

%x = [  x1     x2      x3      x4      x5        x6  ]';
x =  [Pao(1)  Qa(1)  Vve(1)  Pas(1)  Pae(1)   Qvad(1)]';

% Initial states of diodes
Dm = 0; Da = 0;

Da_(1) = Da;
Dm_(1) = Dm;
cae(1) = 0;
rm(1) = 0;

estado_atual = 3;
EDVvec(1) = 0;
ESVvec(1) = 0;
COvec(1) = 0;
EDV = 0;
ESV = 0;
CO = 0;
%% Simulation
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
    end
    if (Dm == 0 && Da ==0) && (estado_atual == 1)
        estado_atual = 2; % se o estado atual é Enchimento(1), vai p/ Contração Isovolumétrica(2)
        %EDV = Vve(i);
    end
    if (Dm == 0 && Da == 1) && (estado_atual == 2)
        estado_atual = 3; % Se o estado atual é Contração Isovolumétrica(2), vai p/ Ejeção(3)
    end
    if (Dm == 0 && Da == 0) && (estado_atual == 3)
        estado_atual = 4; % Se o estado atual é Ejeção(3), vai p/ Contração Isovolumétrica(4)
        %ESV = Vve(i);
        %CO = (EDVvec(i)-ESVvec(i)) * HR
    end
    %COvec(i+1) = CO/1000;
    %EDVvec(i+1) = EDV;
    %ESVvec(i+1) = ESV;
    %SV = EDV - ESV;
    
%     if(estado_anterior ~= estado_atual && estado_anterior == 2)
%         % Encontrar a EDP
%         EDP(i) = 1;
%     else
%         EDP(i) = 0;
%     end

w = (w_rpm(i)/60)*2*pi; % rad/s
% w = (w_rpm/60)*2*pi; % rad/s

% Rk
% if Pve(i) > 1 % x1_
%     Rk = 0;
% else
%     Rk = alfa*(Pve(i) - 1);
% end
% RR = Ri + Ro + Rk + Bo;
RR = Ri + Ro + Bo;

% 0 -> 8000 // Rm = 0.005 e Cae = 30.0

if i < 80000 % 1a constante
    Rm = 0.005;
    Cae = 30.0;
    cae(i+1) = Cae;
    rm(i+1) = Rm;
% 8000 -> 12000 // Rm = Variavel e Cae = variavel
elseif i >= 80000 && i < 120000 % rampa de subida
    Rm = -1e-7*i+0.013000000000000001;
    Cae = 0.01425*i-1110;
    cae(i+1) = Cae;
    rm(i+1) = Rm;
% 12000 -> 32000 // Rm = 0.001 e Cae = 600
elseif i >= 120000 && i < 320000 % 2a constante
    Rm = 0.001;
    Cae = 600;
    cae(i+1) = Cae;
    rm(i+1) = Rm;
% 32000 -> 40000 // Rm = variavel e Cae = 600
elseif i >= 320000 && i <= 400000 % rampa de descida
    Rm = 6.125e-7*i - 0.195;
    Cae = 600;
    cae(i+1) = Cae;
    rm(i+1) = Rm;
end

cae(i+1) = Cae;
rm(i+1) = Rm;
% elseif i >= 40000 % 3a constante
%     Rm = 0.05;
%     Cae = 120.0;
%     cae(i+1) = Cae;
%     rm(i+1) = Rm;
% end

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

x = runkut4(passo,A,x,B);

%w_rpm(i+1) = 100*passo + w_rpm(i); % rpm
phase = 0; % pi/4;
% w_rpm(i+1) = 2000*sin(2*pi*(HR/60)*T(i) + phase) + 10000;
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

PIP(i+1) = PIP(i);

%%
% figure(1)
% plot(T,PIP)
% figure(1)
% subplot(2,1,1)
% plot(T,w_rpm)
% legend('w_{rpm}')
% 
% subplot(2,1,2)
% plot(T,Pao,'k',T,Pve,'r',T,PIP,'b--')
% legend('P_{ao}','P_{ve}','PIP')
% 
% figure(2)
% subplot(1,2,1)
% plot(T,w_rpm)
% 
% subplot(1,2,2)
% plot(T,Qvad)
% 
% figure(2)
% subplot(2,1,1)
% plot(T,Pao,'b',T,Pve,'r',T,Pae,'k')
% legend('Pao','Pve','Pae')
% grid
% 
% subplot(2,1,2)
% plot(T,Qvad,'k')
% legend('Q_{VAD}')
% grid
%%
%plot(T,CO)