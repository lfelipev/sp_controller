clear all
% close all
clc

%% Simulation Time;
start_t = 0;
passo   = 0.0001;
end_t   = 10;

%Uses the already created Time scale
T = start_t:passo:end_t;
n = length(T);

%% Cardiovascular system
HR = 75;
Emax = 2.5;
Emin = 0.06;
En = Elastance(T,passo,HR,end_t);
E = (Emax - Emin)*En + Emin;

% Cardiovascular system model parameters (from Simaan2009);
Rs  = 1.0000; % (0.83-normal,weak; 1.4-severly weak without pump; 0.83-severly weak with pump)(mmHg.sec/mL)
Rm  = 0.0050; % Rm-mitral valve open;(mmHg.sec/mL)
Ra  = 0.0010; % Ra-aortic valve open;(mmHg.sec/mL)
Rc  = 0.0398; % Rc-characteristic resistance;(mmHg.sec/mL)
Cae = 4.4000; % Cr-pulmonary compliance;(mL/mmHg)
Cs  = 1.3300; % Systemic Complinace (ml/mmHg)
Cao = 0.0800; % Aortic Complinace (ml/mmHg)
Ls  = 0.0005; % Ls-inertance of blood in aorta;(mmHg.sec^2/mL)

Vo = 10;

% Preallocating
Pao = zeros(1,length(T));
Qa = zeros(1,length(T));
Vve = zeros(1,length(T));
Pas = zeros(1,length(T));
Pae = zeros(1,length(T));
Pve = zeros(1,length(T));
Dm_ = zeros(1,length(T));
Da_ = zeros(1,length(T));

Rsv = ones(1,length(T));

% Initial Conditions
Pao(1) = 90;
Qa(1)  = 0;
Vve(1) = 140;
Pas(1)  = 90;
Pae(1)  = 5;

Pve(1) = E(1)*(Vve(1) - Vo);

%x = [  x1     x2      x3      x4      x5   ]';
x =  [Pao(1)  Qa(1)  Vve(1)  Pas(1)  Pae(1) ]';

% Initial states of diodes
Dm = 0; Da = 0;

Da_(1) = Da;
Dm_(1) = Dm;

% Início da simulação
for i = 1:n-1
    
% Teste da válvula mitral
if Pae(i) >= Pve(i)
    Dm = 1;
else
    Dm = 0;
end

% Teste da válvula aórtica
if Pve(i) >= Pao(i)
    Da = 1;
else
    Da = 0;
end

% Matrix A
% lambda = E(i) - Emin*En(i);
a13 = (Da/Ra)*(E(i));
a33 = -((Dm/Rm)+(Da/Ra))*E(i);
a53 = (Dm/Rm)*E(i);
a55 = -(1/Rs + Dm/Rm);

%    x = [   x1    x2    x3    x4     x5   ]
%    x = [   Pao   Qa    Vve   Pas    Pae  ]
A(1,:) = [ -Da/Ra  -1    a13    0      0   ]/Cao;
A(2,:) = [    1    -Rc    0    -1      0   ]/Ls;
A(3,:) = [  Da/Ra   0    a33    0    Dm/Rm ];
A(4,:) = [    0     1     0   -1/Rs   1/Rs ]/Cs;
A(5,:) = [    0     0    a53   1/Rs   a55  ]/Cae;

% Matrix B
B(1,1) = (-(Da/Ra)*E(i)*Vo)/Cao;
B(2,1) = 0;
B(3,1) = (Dm/Rm + Da/Ra)*E(i)*Vo;
B(4,1) = 0;
B(5,1) = (-(Dm/Rm)*E(i)*Vo)/Cae;

x = runkut4(passo,A,x,B);

Pao(i+1) = x(1);
Qa(i+1)  = x(2);
Vve(i+1) = x(3);
Pas(i+1) = x(4);
Pae(i+1) = x(5);
Da_(i+1) = Da;
Dm_(i+1) = Dm;

Pve(i+1) = E(i+1)*(Vve(i+1) - Vo);

x1 = 15; y1 = 1.0;
x2 = 25; y2 = 2.0;
m1 = (y2 - y1)/(x2 - x1);

if (i < x1/passo)
    Rs = y1;
elseif (i >= x1/passo) && (i < x2/passo)
    Rs = m1*passo*i - m1*x1 + y1;
else
    Rs = y2;
end

Rsv(i) = Rs;

end

%%
figure(1)
subplot(1,2,1)
plot(T,Pas,'k','LineWidth',2)
ylabel('Pressure (mmHg)')
xlabel('Time (s)')
title('P_s')

subplot(1,2,2)
plot(T,Rsv,'r')
ylabel('Resistance (mmHg.s/mL)')
xlabel('Time (s)')
title('R_s')

