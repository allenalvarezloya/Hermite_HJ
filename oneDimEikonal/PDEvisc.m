function [PDE] = PDEvisc(u,nu,h,m)
ux = (1/h)*[(1:2*m+1)'.*u(2:2*m+2); 0];
uxx = (1/h)^2*[(1:2*m)'.*(2:2*m+1)'.*u(3:2*m+2);0;0];
if ux(1) >= 0
    PDE = -ux + nu*uxx;
else
    PDE = ux + nu*uxx;
end