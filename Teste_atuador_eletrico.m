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
w_rpm = 1000;
w_rpmvec = w_rpm*ones(1, length(t));

J = 0.916e-6;
B = 0.66e-6;
a0 = 0.738e-12;
a1 = 0.198e-10;

Kp = 1.63918032306973e-06;
Ki = 1.85695310261092e-06;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Ew(i) = (w_rpmvec(i)*2*pi/60) - w(i);
    if i ==1
        Eint(i) = 0;
    else
        Eint(i) = Eint(i-1) + Ew(i)*passo;
    end
    
    Te = Kp*Ew(i) + Ki*Eint(i) + (a0*w(i)^3 + a1*Qvad(i)*w(i)^2);

    wdot = xdot_fun_Teste_Te(w(i),Te,Qvad(i));

    w(i+1) = runkut42_Teste_Te(w(i),wdot,Te,Qvad(i),passo);
end

%%
plot(t,w*60/(2*pi), 'LineWidth', 2.5)
grid on
ylabel('\omega (rpm)')
xlabel('time (s)')

% 63% do valor final em t = 0.49s
% t = 1/a
% a = 1/0.49
% tr = 2.2/a = 1.07
% ts = 4/a = 1.96
