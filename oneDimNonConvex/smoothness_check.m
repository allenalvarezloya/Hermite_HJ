n = 80;
xl = -1.0;
xr = 1.0;
m = 2;
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
    lambda = max(abs(sin((1/h)*u(2,:) + 1)));
    dt = lambda*0.5*h;
    if (t + dt > T)
        dt = T - t;
    end
    dth = 0.5*dt;
    ux(1:2*m+1,:) = (1/h)*(1:2*m+1)'.*u(2:2*m+2,:);
    Sd = Sensor(2*m+1,ux,x);
    Vsd = viscosity(Sd);
    v0 = (lambda*h)/(2*m+1);
    for i = 1:n
        u_int = Hmap*[u(1:m+1,i); u(1:m+1,i+1)];
        ud(:,i) = RK4visc(u_int,v0*Vsd(i),dth,h,m);
    end
    ud_ext = [ud(:,end) ud ud(:,1)];
    udx(1:2*m+1,:) = 1/h*(1:2*m+1)'.*ud_ext(2:2*m+2,:);
    Sp = Sensor(2*m+1,udx,xd_ext);
    Vsp = viscosity(Sp);
    for i = 1:n+1
        u_int = Hmap*[ud_ext(1:m+1,i); ud_ext(1:m+1,i+1)];
        u(:,i) = RK4visc(u_int,v0*Vsp(i),dth,h,m);
    end
    t = t + dt;
    TITLE = sprintf('time = %3.2f',t);
    subplot(2,1,1)
    plot(x,u(1,:),xd,Vsd,'LineWidth',2);
    title(TITLE)
    axis([-1.0 1.0 -1.1 1.5])
    subplot(2,1,2)
    plot(xd,Sd,'LineWidth',2)
    title(TITLE)
    axis([-1.0 1.0 0.0 13])
    pause(0.5)
end

% Refined solution 
nr = 4; % Number of grid refinement has to be even
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
h_refined = (xr - xl)/N;
x_refined = xl + (0:N)*h_refined;
u_actual = zeros(N+1,1);



P = @(z) pi*sin(pi*z);

U_exact = zeros(length(x_refined),1);
prev = zeros(length(x_refined),1);

for i = 1:length(x_refined)
%     x_coord = Newton_solve(x_refined(i),prev(i),T);
%     prev(i) = x_coord;
    f = @(z) t*sin(pi*sin(pi*z)+1) + z - x_refined(i);
    x_coord = fzero(f,0);
    f(x_coord)
    p = P(x_coord);
    u_actual(i) = T*(sin(p+1)*p + cos(p+1)) - cos(pi*x_coord); 
end
plot(x_refined,u_refined,'LineWidth',2)
hold on
plot(x_refined(1:12:end),u_actual(1:12:end),'o')
hold off
legend('Numerical Solution','Actual Solution')
saveas(gcf,'Example_03_1D','epsc')


function [x_sol] =  Newton_solve(x,x_guess,t)
    f = @(z) t*sin(pi*sin(pi*z)+1) + z - x;
    fp = @(z) t*cos(pi*sin(pi*z)+1)*pi^2*cos(pi*z) + 1.0;

    tol = 10e-16;
    while abs(f(x_guess)) > tol
        dx = -f(x_guess)/fp(x_guess);
        x_guess = x_guess + dx;
        abs(f(x_guess))
    end
    x_sol = x_guess;
end


