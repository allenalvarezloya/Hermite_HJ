% This script computes and plots the numerical solution
% of the Eikonal equation.
clear, clc;
n = 41;
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
    Sd = Sensor(2*m+1,m,u,x);
    Vsd = viscosity(Sd);
    lambda = 1.0;
    v0 = (lambda*h)/(2*m+1);
    for i = 1:n
        u_int = Hmap*[u(1:m+1,i); u(1:m+1,i+1)];
        ud(:,i) = RK4visc(u_int,v0,dth,h,m);
    end
    ud_ext = [ud(:,end) ud ud(:,1)];
    Sp = Sensor(2*m+1,m,ud_ext,xd_ext);
    Vsp = viscosity(Sp);
    for i = 1:n+1
        u_int = Hmap*[ud_ext(1:m+1,i); ud_ext(1:m+1,i+1)];
        u(:,i) = RK4visc(u_int,v0,dth,h,m);
    end
    t = t + dt;
    % Compute numerical solution on refined grid  
    nr = 8; % Number of grid refinement has to be even
    grid_refine = (-0.5:1/nr:0.5)';
    Eval_Mat = zeros(nr,2*m+2);
    Eval_Mat_Left = zeros(nr/2+1,2*m+2);
    Eval_Mat_Right = zeros(nr/2,2*m+2);
    for idx = 0:2*m+1
        Eval_Mat(:,1+idx) = grid_refine(2:end).^(idx);
        Eval_Mat_Left(:,1+idx) = grid_refine(nr/2+1:end).^(idx);
        Eval_Mat_Right(:,1+idx) = grid_refine(2:nr/2+1).^(idx);
    end
    u_refined = [];
    u_refined = Eval_Mat_Left*u(:,1);
    for k = 2:n
        u_refined = [u_refined; Eval_Mat*u(:,k)];
    end
    u_refined = [u_refined; Eval_Mat_Right*u(:,n+1)];
    N = length(u_refined)-1;
    h_refined = 2*pi/N;
    x_refined = 0 + (0:N)*h_refined;
    TITLE = sprintf('time = %3.2f',t);
    plot(x_refined,u_refined,'LineWidth',2);
    axis([0 2*pi -1 1])
    legend('Numerical Solution')
    title(TITLE)
    drawnow
end


