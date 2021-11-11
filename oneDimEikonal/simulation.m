% This script computes and plots the numerical solution
% of the Eikonal equation
clear, clc;
n = 80;
xl = 0.0;
xr = 2*pi;
m = 3;
h = (xr - xl)/n;
x = xl + h*(0:n)';
xd = xl + h*(1:n)' - 0.5*h;
xd_ext = [xd(end); xd; xd(1)];
ud = zeros(2*m+2,n);
u = initialData(x,h,m);
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
    Sd = Sensor(2*m+3,m,u,x);
    Vsd = viscosity(Sd);
    lambda = 1.0;
    v0 = (lambda*h)/(2*m+1);
    for i = 1:n
        u_int = Hmap*[u(1:m+1,i); u(1:m+1,i+1)];
        ud(:,i) = RK4visc(u_int,v0,dth,h,m);
    end
    ud_ext = [ud(:,end) ud ud(:,1)];
    Sp = Sensor(2*m+3,m,ud_ext,xd_ext);
    Vsp = viscosity(Sp);
    for i = 1:n+1
        u_int = Hmap*[ud_ext(1:m+1,i); ud_ext(1:m+1,i+1)];
        u(:,i) = RK4visc(u_int,v0,dth,h,m);
    end
    t = t + dt;
    TITLE = sprintf('time = %3.2f',t);
    plot(x,u(1,:),'LineWidth',2);
    axis([0 2*pi -1 1])
    legend('Numerical Solution')
    title(TITLE)
    drawnow
end


