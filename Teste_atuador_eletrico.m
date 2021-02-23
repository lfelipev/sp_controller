clear
% close all
clc

global B a0 a1 J

%% Simulation Time;
start_t = 0;
passo   = 0.0001;
end_t   = 10;

% Uses the already created Time scale
t = start_t:passo:end_t;
n = length(t);

Qvad = 100 + 92*sin(2*pi*(1/0.67)*t);

%  constant speed
w    = zeros(1, length(t));
Eint = zeros(1, length(t));
Ew   = zeros(1, length(t));

J = 0.916e-6;
B = 0.66e-6;
a0 = 0.738e-12;
a1 = 0.198e-10;

Kp = 2*B*10;
Ki = Kp*B/J;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Ew(i) = 100 - w(i);
    if i ==1
        Eint(i) = 0;
    else
        Eint(i) = Eint(i-1) + Ew(i)*passo;
    end
    Te = Kp*Ew(i) + Ki*Eint(i);

    wdot = xdot_fun_Teste_Te(w(i),Te,Qvad(i));

    w(i+1) = runkut42_Teste_Te(w(i),wdot,Te,Qvad(i),passo);
end

plot(t,w)
