
function x = runkut4(passo,A,x,B)

xdot = A*x + B;
kx1 = passo*xdot;

x1 = x + 0.5*kx1;
xdot = A*x1 + B;
kx2 = passo*xdot;

x1 = x + 0.5*kx2;
xdot = A*x1 + B;
kx3 = passo*xdot;

x1 = x + kx3;
xdot = A*x1 + B;
kx4 = passo*xdot;

x = x + (kx1 + 2*kx2 + 2*kx3 + kx4)/6;