clear all
% close all
clc

%% Simulation Time;
start_t = 0;
passo   = 0.0001;
end_t   = 5;
HR = 75;

t = start_t:passo:end_t;
n = length(t);

for i = 1:n
    w_rpm(i) = 1000*sin(2*pi*(HR/60)*t(i)) + 9000;
end

plot(t,w_rpm)