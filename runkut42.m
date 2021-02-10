function x = runkut42(x,xdot,E,w,passo)

kx1 = passo*xdot;

x1 = x + 0.5*kx1;

xdot = xdot_fun(x1,E,w);
kx2 = passo*xdot;

x1 = x + 0.5*kx2;
xdot = xdot_fun(x1,E,w);
kx3 = passo*xdot;

x1 = x + kx3;
xdot = xdot_fun(x1,E,w);
kx4 = passo*xdot;

x = x + (kx1 + 2*kx2 + 2*kx3 + kx4)/6;