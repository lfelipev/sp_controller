function xdot = xdot_fun_Te(x,E,Te)

global Rs Ra Rm Rc Cao Cs Cae Ls Dm Da Vo RR LL B2 B a0 a1 J

% Matrix A
% lambda = E(i) - Emin*En(i);
a13 = (Da/Ra)*E;
a33 = -((Dm/Rm)+(Da/Ra))*E;
a53 = (Dm/Rm)*E;
a55 = -(1/Rs + Dm/Rm);

xdot(1,1) = (-Da/Ra*x(1)   -1*x(2)   +a13*x(3)    +0*x(4)      +0*x(5)    +1*x(6) -    (Da/Ra)*E*Vo     )/Cao;
xdot(2,1) = (     1*x(1)  -Rc*x(2)     +0*x(3)    -1*x(4)      +0*x(5)    +0*x(6) +           0            )/Ls;
xdot(3,1) = ( Da/Ra*x(1)   +0*x(2)   +a33*x(3)    +0*x(4)  +Dm/Rm*x(5)    -1*x(6) + (Dm/Rm + Da/Ra)*E*Vo);
xdot(4,1) = (     0*x(1)   +1*x(2)     +0*x(3) -1/Rs*x(4)   +1/Rs*x(5)    +0*x(6) +           0             )/Cs;
xdot(5,1) = (     0*x(1)   +0*x(2)   +a53*x(3) +1/Rs*x(4)    +a55*x(5)    +0*x(6) -    (Dm/Rm)*E*Vo      )/Cae;
xdot(6,1) = (    -1*x(1)   +0*x(2)     +E*x(3)    +0*x(4)      +0*x(5)   -RR*x(6) -   (E*Vo) - B2*x(7)^2    )/LL;
xdot(7,1) = (Te - B*x(7) - (a0*x(7)^3 + a1*x(6)*x(7)^2))/J;




