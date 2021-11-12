% This script computes and plots the numerical
% Solution of a one-dimensional Riemann problem
clear, clc;
n = 80;
xl = -1.0;
xr = 1.0;
m = 3;
h = (xr - xl)/n;
x = xl + h*(0:n)';
xd = xl + h*(1:n)' - 0.5*h;
xd_ext = [xd(end); xd; xd(1)];
ud = zeros(2*m+2,n);
u1 = initialData(x,h,m,1);
u2 = initialData(x,h,m,2);
u = 0.0*u1;
ux = zeros(2*m+2,n+1);
udx = zeros(2*m+2,n+2);
Hmap = Hermite_map(m,0.0,1.0,0.5,0);
t = 0.0;
T = 1.0;
while t < T
    dt = 0.1*h;
    if (t + dt > T)
        dt = T - t;
    end
    dth = 0.5*dt;
    Sd = Sensor(2*m+1,u,x);
    Vsd = viscosity(Sd);
    lambda = 2.0;
    v0 = (lambda*h)/(2*m+1);
    for i = 1:0.5*n
        u_int = Hmap*[u1(1:m+1,i); u1(1:m+1,i+1)];
        ud(:,i) = RK4visc(u_int,v0*Vsd(i),dth,h,m);
    end
    for i = 0.5*n+1:n
        u_int = Hmap*[u2(1:m+1,i); u2(1:m+1,i+1)];
        ud(:,i) = RK4visc(u_int,v0*Vsd(i),dth,h,m);
    end
    ud_ext = [ud(:,end) ud ud(:,1)];
    Sp = Sensor(2*m+1,ud,xd);
    Vsp = viscosity(Sp);
    Vsp = [0.0; Vsp; 0.0];
    for i = 2:n
        u_int = Hmap*[ud_ext(1:m+1,i); ud_ext(1:m+1,i+1)];
        u(:,i) = RK4visc(u_int,v0*Vsp(i),dth,h,m);
    end
    u(:,1) = [-2.0; 2*h; zeros(2*m,1)]; 
    u(:,end) = [-2.0; -2*h; zeros(2*m,1)];
    u1 = u;
    u2 = u;
    t = t + dt;
    TITLE = sprintf('time = %3.2f',t);
    plot(x,u(1,:),'LineWidth',2);
    title(TITLE)
    axis([-1.0 1.0 -2.0 0.0])
    drawnow
end
