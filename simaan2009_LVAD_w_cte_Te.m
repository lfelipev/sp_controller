clear
% close all
clc

global Da Dm Rs Rm Ra Rc Cae Cs Cao Ls RR LL Vo B2 alpha B a0 a1 J

%% Simulation Time;
start_t = 0;
passo   = 0.0001;
end_t   = 30;

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

% Pump parameters
Bo = 0.296;
B1 = 0.027;
B2 = -9.33e-5;

LL = Li + Lo + B1;
alpha = -3.5;

Vo = 10;

Pao = zeros(1, length(T));
Qa  = zeros(1, length(T));
Vve = zeros(1, length(T));
Pas = zeros(1, length(T));
Pae = zeros(1, length(T));
Qvad = zeros(1, length(T));

%  constant speed
w_rpm = 10000;
wrads = (w_rpm*2*pi/60);
w = wrads*ones(1, length(T));

Pve  = zeros(1, length(T));
Eint = zeros(1, length(T));
Ew = zeros(1, length(T));

% Initial Conditions
Pao(1)  = 80;
Qa(1)   = 0;
Vve(1)  = 100;
Pas(1)  = 75;
Pae(1)  = 16;
Qvad(1) = 40; % x6 - LVAD flow

%x = [  x1     x2      x3      x4      x5        x6      x7  ]';
x =  [Pao(1)  Qa(1)  Vve(1)  Pas(1)  Pae(1)   Qvad(1)   w(1) ]';

% Initial states of diodes
Dm = 0; Da = 0;

Tevec = zeros(1, length(T));

J = 0.916e-6;
B = 0.66e-6;
a0 = 0.738e-12;
a1 = 0.198e-10;

% Ganhos sintonizados pelo MATLAB
Kp = 1.63918032306973e-06;
Ki = 1.85695310261092e-06;

% Ganhos para cancelamento de Pólos 
% Kp = 8*B;
% Ki = Kp*B/J;

% Contador percentual da simula��o
ip = 0;
lg = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lógica dos diodos
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

% Resistor do LVAD
RR = Ri + Ro + Bo;

% Variação da referência em 15s
if i >= 150000
    w_rpm = 13000;
end

% Calculo do erro e da integral do erro
Ew(i) = (w_rpm*2*pi/60) - w(i);
if i ==1
    Eint(i) = 0;
else
    Eint(i) = Eint(i-1) + Ew(i)*passo;
end

% Cálculo do torque elétrico
% Te = Kp*Ew(i);
Te = Kp*Ew(i) + Ki*Eint(i) + (a0*w(i)^3 + a1*Qvad(i)*w(i)^2);
Tevec(i) = Te;

% Função das variáveis de estado
xdot = xdot_fun_Te(x,E(i),Te);

% Runge-Kutta
x = runkut42_Te(x,xdot,E(i),Te,passo);

Pao(i+1) = x(1);
Qa(i+1)  = x(2);
Vve(i+1) = x(3);
Pas(i+1) = x(4);
Pae(i+1) = x(5);
Qvad(i+1)= x(6);
w(i+1)   = x(7);

Pve(i+1) = E(i+1)*(Vve(i+1) - Vo);


if (i > ip)
    fprintf('Executing ... \t%d %%\r',lg);
    lg = lg + 10;
    ip = ip + (n-1)/10;
end

end
%%
figure(1)
subplot(3,1,1)
plot(T,Qvad)
grid on
title('Fluxo através do VAD (ml/s)')
subplot(3,1,2)
plot(T,w)
grid on
title('Velocidade de rotação do VAD (rad/s)')
subplot(3,1,3)
plot(T,Tevec)
title('Torque elétrico (N/m)')
grid on