function [PIP, SP, COvec, EDP, Pas, rs] = physiological_simaan(enable_preload, end_t)

global Rs Ra Rm Rc Cao Cs Cae Ls Dm Da Vo

% clear all
% close all
clc

%% Simulation Time;
start_t = 0;
passo   = 0.0001;

%Uses the already created Time scale
T = start_t:passo:end_t;
n = length(T);

%% Cardiovascular system
HR = 75;
Emax = 2.00;
Emin = 0.05;
En = Elastance(T,passo,HR,end_t);
E = (Emax - Emin)*En + Emin;

% Cardiovascular system model parameters (from Simaan2009);
Rs  = 1.0000; % (0.83-normal,weak; 1.4-severly weak without pump; 0.83-severly weak with pump)(mmHg.sec/mL)
Rm  = 0.0050; % Rm-mitral valve open;(mmHg.sec/mL)
Cae = 4.4000; % Cr-pulmonary compliance;(mL/mmHg)
Ra  = 0.0010; % Ra-aortic valve open;(mmHg.sec/mL)
Rc  = 0.0398; % Rc-characteristic resistance;(mmHg.sec/mL)
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
rs = ones(1, length(T));

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

CO = 0;
EDV = 0;
ESV = 0;
estado_atual = 3;
d_PIP(1) = 0;
SP      = zeros(1,length(T));
EDP(1) = Pve(1);
Aux(1) = 40;
estado(1) = 0;

% Início da simulação
w = 0; % não existe velocidade aqui!
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

% Detector de estados do ciclo cardiaco
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

xdot = xdot_fun_phy(x,E(i)); % não existe velocidade aqui!
x = runkut42_phy(x,xdot,E(i),passo);

PIP(i) = Pve(i);

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
    Aux(i) = 40;
else
    EDP(i) = EDP(i-1);
    Aux(i) = 30;
end
end

Pao(i+1) = x(1);
Qa(i+1)  = x(2);
Vve(i+1) = x(3);
Pas(i+1) = x(4);
Pae(i+1) = x(5);

Pve(i+1) = E(i+1)*(Vve(i+1) - Vo);
end
PIP(i+1) = 0;
EDP(i+1) = 0;
end