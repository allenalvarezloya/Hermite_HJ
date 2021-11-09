function [out] = PDEvisc(u,nu,h,m)
ux = (1/h)*[(1:2*m+1)'.*u(2:2*m+2); 0];
uxx = (1/h)^2*[(1:2*m)'.*(2:2*m+1)'.*u(3:2*m+2);0;0];
ONE = [1; zeros(2*m+1,1)];
FOUR = [4; zeros(2*m+1,1)];
uxux = conv(ux,ux);
RHS = -0.25*conv(uxux(1:2*m+2) - ONE,uxux(1:2*m+2) - FOUR);
out = RHS(1:2*m+2)+nu*uxx;

end