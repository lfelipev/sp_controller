clear
close all
clc

global Rs Ra Rm Rc Cao Cs Cae Ls Dm Da Vo

%% Simulation Time;
start_t = 0;
passo   = 0.0001;
end_t   = 10;

%Uses the already created Time scale
T = start_t:passo:end_t;
n = length(T);

%% Cardiovascular system
HR = 75;
Emax = 2;
Emin = 0.05;
En = Elastance(T,passo,HR,end_t);
E = (Emax - Emin)*En + Emin;

% Cardiovascular system model parameters (from Simaan2009);
Rs  = 0.500; % (0.83-normal,weak; 1.4-severly weak without pump; 0.83-severly weak with pump)(mmHg.sec/mL)
Rm  = 0.1; % Rm-mitral valve open;(mmHg.sec/mL)
Ra  = 0.0010; % Ra-aortic valve open;(mmHg.sec/mL)
Rc  = 0.0398; % Rc-characteristic resistance;(mmHg.sec/mL)
Cae = 4.4000; % Cr-pulmonary compliance;(mL/mmHg)
Cs  = 1.3300; % Systemic Complinace (ml/mmHg)
Cao = 0.0800; % Aortic Complinace (ml/mmHg)
Ls  = 0.0005; % Ls-inertance of blood in aorta;(mmHg.sec^2/mL)

Vo = 10;

% Preallocating
SP     = zeros(1,length(T));
rs     = ones(1, length(T));
rm     = zeros(1, length(T));
cae    = zeros(1, length(T));
Pao    = zeros(1, length(T));
Qa     = zeros(1, length(T));
Vve    = zeros(1, length(T));
Pas    = zeros(1, length(T));
Pae    = zeros(1, length(T));
Qvad   = zeros(1, length(T));
Pve    = zeros(1, length(T));
PIP    = zeros(1, length(T));
d_PIP  = zeros(1, length(T));
estado = zeros(1, length(T));
COvec_phy  = zeros(1, length(T));

% Initial Conditions
Pao(1) = 80;
Qa(1)  = 0;
Vve(1) = 100;
Pas(1)  = 75;
Pae(1)  = 16;

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
PIP(1) = Pve(1);
estado(1) = 0;
enable_preload = 1;
SV = 0;

% Contador percentual da simulação
ip = 0;
lg = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n-1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
        SV = EDV-ESV;
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
        esv(i) = Vve(i);
    end
    
    COvec_phy(i+1) = ((SV) * HR) / 1000;
    %COvec_phy(i+1) = Qa(i);
    
    estado(i+1) = estado_atual;
    
    a13 = (Da/Ra)*E(i);
    a33 = -((Dm/Rm)+(Da/Ra))*E(i);
    a53 = (Dm/Rm)*E(i);
    a55 = -(1/Rs + Dm/Rm);
    
    xdot = xdot_fun_phy(x,E(i));
    x = runkut42_phy(x,xdot,E(i),passo);
    
    if i > 1
        d_PIP(i) = (PIP(i) - PIP(i-1))/passo;
        if d_PIP(i-1) >= 0.01 && d_PIP(i) < -0.01
            SP(i) = PIP(i);
        else
            SP(i) = SP(i-1);
        end
    end
    
    Pao(i+1) = x(1);
    Qa(i+1)  = x(2);
    Vve(i+1) = x(3);
    Pas(i+1) = x(4);
    Pae(i+1) = x(5);
    
    Pve(i+1) = E(i+1)*(Vve(i+1) - Vo);
    PIP(i+1) = Pve(i+1);
    
    if (i > ip)
        fprintf('Executing ... \t%d %%\r',lg);
        lg = lg + 10;
        ip = ip + (n-1)/10;
    end
end
%%
Vve_05 = Vve;
Pve_05 = Pve;
index_05 = length(esv);
save('preload_variation_05.mat','Vve_05','Pve_05', 'index_05')
