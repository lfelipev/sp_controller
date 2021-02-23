function xdot = xdot_fun_Teste_Te(w,Te,Qvad)

global B a0 a1 J

xdot = (Te - B*w - (a0*w^3 + a1*Qvad*w^2))/J;




