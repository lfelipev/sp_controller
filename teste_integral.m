clear
close all
clc

passo = 0.001;
t = 0:passo:1;
N = length(t);
x = sin(2*pi*t);
% x = ones(1,length(t));

plot(t,x)
hold on

intx = zeros(1,length(t));
for i = 1:N-1
    intx(i+1) = intx(i) + x(i)*passo;
end
plot(t,intx)