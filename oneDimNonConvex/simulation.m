% This matlab script computes and plots the 
% numerical solution of a on dimensional 
% Hamilton-Jacobi problem with a nonconvex Hamiltonian

n = 80;
xl = -1.0;
xr = 1.0;
m = 3;
h = (xr - xl)/n;
x = xl + h*(0:n)';
xd = xl + h*(1:n)' - 0.5*h;
xd_ext = [xd(end); xd; xd(1)];
ud = zeros(2*m+2,n);
u = initialData(x,h,m);
ux = zeros(2*m+2,n+1);
udx = zeros(2*m+2,n+2);
Hmap = Hermite_map(m,0.0,1.0,0.5,0);
t = 0.0;
T = 1.5/pi^2;
while t < T
    dt = 0.1*h;
    if (t + dt > T)
        dt = T - t;
    end
    dth = 0.5*dt;
    Sd = Sensor(2*m+3,m,u,x);
    Vsd = viscosity(Sd);
    lambda = max(abs(sin((1/h)*u(2,:) + 1)));
    v0 = (lambda*h)/(2*m+1);
    for i = 1:n
        u_int = Hmap*[u(1:m+1,i); u(1:m+1,i+1)];
        ud(:,i) = RK4visc(u_int,v0*Vsd(i),dth,h,m);
    end
    ud_ext = [ud(:,end) ud ud(:,1)];
    Sp = Sensor(2*m+3,m,ud_ext,xd_ext);
    Vsp = viscosity(Sp);
    for i = 1:n+1
        u_int = Hmap*[ud_ext(1:m+1,i); ud_ext(1:m+1,i+1)];
        u(:,i) = RK4visc(u_int,v0*Vsp(i),dth,h,m);
    end
    t = t + dt;
    TITLE = sprintf('time = %3.2f',t);
    plot(x,u(1,:),'LineWidth',2);
    title(TITLE)
    drawnow
    t = t + dt;
end
