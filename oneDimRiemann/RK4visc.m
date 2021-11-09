function [u_p] = RK4visc(u,nu,dt,h,m)
K1 = PDEvisc(u,nu,h,m);
K2 = PDEvisc(u + 0.5*dt*K1,nu,h,m);
K3 = PDEvisc(u + 0.5*dt*K2,nu,h,m);
K4 = PDEvisc(u + dt*K3,nu,h,m);

u_p = u + (dt/6)*(K1 + 2*K2 + 2*K3 + K4);
end