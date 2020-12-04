function [PIP, SP, COvec, EDP] = physiological_simaan(enable_preload, end_t)

    %clear all
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
    %E(26000:36000) = (1)*En(30000:40000) + Emin; PVC

    % Cardiovascular system model parameters (from Simaan2009);

    % Rs  = 1.0000; % (0.83-normal,weak; 1.4-severly weak without pump; 0.83-severly weak with pump)(mmHg.sec/mL)
    % Rm  = 0.0050; % Rm-mitral valve open;(mmHg.sec/mL)

    Rs = 1;
    Rm = 0.1;

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
            Rm = -2.375e-6*i+0.48;
            %Cae = 1.8e-3*i -136;
            cae(i+1) = Cae;
            rm(i+1) = Rm;
        % 12000 -> 32000 // Rm = 0.001 e Cae = 600
        elseif i >= 200000 && i < 400000 % 2a constante
            Rm = 0.005;
            %Cae = 80;
            cae(i+1) = Cae;
            rm(i+1) = Rm;
        % 32000 -> 40000 // Rm = variavel e Cae = 600
        elseif i >= 400000 && i <= 500000 % rampa de descida
            Rm = 9.5e-7*i - 0.375;
            %Cae = 80;
            cae(i+1) = Cae;
            rm(i+1) = Rm;
        elseif i >= 500000 && i <= 700000
            Rm = 0.1;
            %Cae = 80;
            cae(i+1) = Cae;
            rm(i+1) = Rm;
        elseif i >= 700000 && i <= 800000
            Rm = 1.5e-6*i - 0.95;
            rm(i+1) = Rm;
        elseif i>=800000
            Rm = 0.25;
            rm(i+1) = Rm;
        end
    end
%         % Calculo da Preload
%         if enable_preload
%             % Varia?ão da Pré-carga
%             if i < 80000 % 1a constante
%                 Rm = 0.1;
%                 cae(i+1) = Cae;
%                 rm(i+1) = Rm;
%             % 8000 -> 12000 // Rm = Variavel e Cae = variavel
%             elseif i >= 80000 && i < 120000 % rampa de subida
%                 Rm = -2.375e-6*i+0.29000000000000004;
%                 %Cae = 1.8e-3*i -136;
%                 cae(i+1) = Cae;
%                 rm(i+1) = Rm;
%             % 12000 -> 32000 // Rm = 0.001 e Cae = 600
%             elseif i >= 120000 && i < 320000 % 2a constante
%                 Rm = 0.005;
%                 %Cae = 80;
%                 cae(i+1) = Cae;
%                 rm(i+1) = Rm;
%             % 32000 -> 40000 // Rm = variavel e Cae = 600
%             elseif i >= 320000 && i <= 400000 % rampa de descida
%                 Rm = 3.0625e-6*i - 0.9749999999999999;
%                 %Cae = 80;
%                 cae(i+1) = Cae;
%                 rm(i+1) = Rm;
%             elseif i >= 400000
%                 Rm = 0.25;
%                 %Cae = 80;
%                 cae(i+1) = Cae;
%                 rm(i+1) = Rm;
%             end
%         end

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

        % Calculo do PIP que nesse caso é o Pve

        PIP(i) = Pve(i);
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

        Pao(i+1) = x(1);
        Qa(i+1)  = x(2);
        Vve(i+1) = x(3);
        Pas(i+1) = x(4);
        Pae(i+1) = x(5);
        Da_(i+1) = Da;
        Dm_(i+1) = Dm;

        Pve(i+1) = E(i+1)*(Vve(i+1) - Vo);
    end
    PIP(i+1) = 0;
    EDP(i+1) = 0;
end