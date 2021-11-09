function [PDE] = PDEvisc(u,nu,h,m)
ux = (1/h)*[(1:2*m+1)'.*u(2:2*m+2); 0];
uxx = (1/h)^2*[(1:2*m)'.*(2:2*m+1)'.*u(3:2*m+2);0;0];
uxux = conv(ux,ux);
PDE = -0.5*uxux(1:2*m+2) + nu*uxx;
end