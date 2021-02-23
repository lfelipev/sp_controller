function x = runkut42_Teste_Te(x,xdot,Te,Qvad,passo)

kx1 = passo*xdot;

x1 = x + 0.5*kx1;

xdot = xdot_fun_Teste_Te(x,Te,Qvad);
kx2 = passo*xdot;

x1 = x + 0.5*kx2;
xdot = xdot_fun_Teste_Te(x,Te,Qvad);
kx3 = passo*xdot;

x1 = x + kx3;
xdot = xdot_fun_Teste_Te(x,Te,Qvad);
kx4 = passo*xdot;

x = x + (kx1 + 2*kx2 + 2*kx3 + kx4)/6;