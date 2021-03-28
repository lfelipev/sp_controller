function RMSE = objectiveFunction(delta)
% clear
% close all
% clc

global Rs Ra Rm Rc Cao Cs Cae Ls Dm Da Vo RR LL B2 B a0 a1 J
tic
%% Simulation Time;
start_t = 0;
passo   = 0.0001;
end_t   = 120;

%Uses the already created Time scale
T = start_t:passo:end_t;
n = length(T);

%% Cardiovascular system
HR = 100;
Emax = 1.2;
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

% Pump parameters
Bo = 0.17070;
B1 = 0.02177;
B2 = -9.9025e-5;

LL = Li + Lo + B1;
alpha = -3.5;

Vo = 10;

SP  = zeros(1,length(T));
PIP = zeros(1, length(T));
d_PIP  = zeros(1, length(T));
COvec  = zeros(1, length(T));
EDVvec = zeros(1, length(T));
ESVvec = zeros(1, length(T));
rs  = ones(1, length(T));
cae = zeros(1, length(T));
Pao = zeros(1, length(T));
Qa  = zeros(1, length(T));
Vve = zeros(1, length(T));
Pas = zeros(1, length(T));
Pae = zeros(1, length(T));
Qvad = zeros(1, length(T));
Pve  = zeros(1, length(T));
Ew = zeros(1, length(T));
Eint = zeros(1, length(T));
Tevec = zeros(1, length(T));
estado = zeros(1,length(T));
Qao = zeros(1,length(T));
Dmvec = zeros(1,length(T));
Davec = zeros(1,length(T));
Cycle = zeros(1,length(T));
kspvec = zeros(1,length(T));

%  constant speed
w_rpm = 11500;
w_rpmvec = w_rpm*ones(1, length(T));
wrads = (w_rpm*2*pi/60);
w = wrads*ones(1, length(T));

% Initial Conditions
Pao(1)  = 80;
Qa(1)   = 0;
Vve(1)  = 100;
Pas(1)  = 75;
Pae(1)  = 16;
Qvad(1) = 40; % x6 - LVAD flow
Qao(1) = 0;
Dmvec(1) = 0;
Davec(1) = 0;
CO = 0;
EDV = 0;
ESV = 0;
Pve(1) = E(1)*(Vve(1) - Vo);
kspvec(1) = 0;
estado_atual = 3;

%x = [  x1     x2      x3      x4      x5        x6      x7 ]';
x =  [Pao(1)  Qa(1)  Vve(1)  Pas(1)  Pae(1)   Qvad(1)   w(1)]';

% Initial states of diodes
Dm = 0; Da = 0;

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

% Contador percentual da simulação
ip = 0;
lg = 0;

Nref = 11500;
SPref = 83.55;
ksp = 10;
COvec_aux = 0;

Qvad_int = 0;
count = 1;

counter_erro = 1;
%0.1 segundos
counter_error_constant = 10000;
COerr = 0;

load('simaan2009_phy.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n-1
    counter_erro = counter_erro+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if i > count*60/HR/passo
        count = count + 1;
        Cycle(i) = 10;
    end
    
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
%     
%     j = j + 1;
%     if j >= round(0.6667/passo)
%         COvec_aux = mean(Qvad(i-j+1:i));
%         j = 0;
%     end
%     COvec(i) = COvec_aux * 0.06;
    
    % Varia?ão da Pré-carga
    if i < 160000 % 1a constante
        Rm = 0.1;
        cae(i+1) = Cae;
        % 8000 -> 12000 // Rm = Variavel e Cae = variavel
    elseif i >= 160000 && i < 200000 % rampa de subida
        Rm = -5e-7*i+0.18;
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
    
    %Rk
    if Pve(i) > 1 % x1_
        Rk = 0;
    else
        Rk = 0;
        %Rk = alpha*(Pve(i) - 1);
    end
    
    RR = Ri + Ro + Rk + Bo;
    %RR = Ri + Ro + Bo;
    
    % Variação da referência em 15s
    %     if i >= 150000
    %         w_rpm = 13000;
    %     end
    
    % Calculo do erro e da integral do erro
    Ew(i) = (w_rpmvec(i)*2*pi/60) - w(i);
    if i ==1
        Eint(i) = 0;
    else
        Eint(i) = Eint(i-1) + Ew(i)*passo;
    end
    
    % Cálculo do torque elétrico
    Te = Kp*Ew(i) + Ki*Eint(i) + (a0*w(i)^3 + a1*Qvad(i)*w(i)^2);
    %Te = Kp*Ew(i) + Ki*Eint(i);
    Tevec(i) = Te;
    
    % Função das variáveis de estado
    xdot = xdot_fun_Te(x,E(i),Te);
    
    PIP(i) = Pve(i) - Li*xdot(6,1) - Ri*x(6);
    if i > 1
        d_PIP(i) = (PIP(i) - PIP(i-1))/passo;
        if d_PIP(i-1) >= 0.01 && d_PIP(i) < -0.01
            SP(i) = PIP(i);
        else
            SP(i) = SP(i-1);
        end
    end
%     

    if(counter_erro == counter_error_constant)
        counter_erro = 1;
        COerr = COvec_phy(i-1) - COvec(i-1);
        
        if COerr >= 0
            ksp = ksp - delta;
        elseif COerr < 0 && COerr > 0.2
            ksp = ksp + delta;
        else
            ksp = ksp + 1;
        end
    end
    
    % limite inferior de ksp
    if ksp < 5
        ksp = 5;
    end
    % limite superior de ksp
    if ksp > 150
        ksp = 150;
    end
    
    kspvec(i+1) = ksp;
    w_rpmvec(i+1) = ksp*(SP(i)-SPref) + Nref;
    
    % Runge-Kutta
    x = runkut42_Te(x,xdot,E(i),Te,passo);
    
    Pao(i+1) = x(1);
    Qa(i+1)  = x(2);
    Vve(i+1) = x(3);
    Pas(i+1) = x(4);
    Pae(i+1) = x(5);
    Qvad(i+1)= x(6);
    w(i+1) = x(7);
    
    Pve(i+1) = E(i+1)*(Vve(i+1) - Vo);
    Qao(i+1) = (Pao(i+1) - Pve(i+1))/Ra;
    
    if (i > ip)
        fprintf('Executing ... \t%d %%\r',lg);
        lg = lg + 10;
        ip = ip + (n-1)/10;
    end
    
    Qvad_int = Qvad_int + Qvad(i)*passo;
    
    if Cycle(i) == 10
        CO = Qvad_int;
        Qvad_int = Qvad_int*0;
    end
    COvec(i) = CO*0.06;
end
PIP(i+1) = PIP(i);
kspvec(i+1) = ksp;
toc

COvec_vg = COvec;
SP_vg = SP;
save('simaan2009_LVAD_VG_Cont.mat','COvec_vg','SP_vg')
load('simaan2009_con.mat')
load('simaan2009_LVAD_SP_Cont.mat')
load('simaan2009_phy.mat')

RMSE = sqrt(mean((COvec_vg-COvec_phy).^2));

