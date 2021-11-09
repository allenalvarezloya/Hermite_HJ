function [u_p] = RK4visc(u,nu,dt,h,m,Cx,Sx,Sy)
K1 = PDEvisc(u,nu,h,m,Cx,Sx,Sy);
K2 = PDEvisc(u + 0.5*dt*K1,nu,h,m,Cx,Sx,Sy);
K3 = PDEvisc(u + 0.5*dt*K2,nu,h,m,Cx,Sx,Sy);
K4 = PDEvisc(u + dt*K3,nu,h,m,Cx,Sx,Sy);

u_p = u + (dt/6)*(K1 + 2*K2 + 2*K3 + K4);
end