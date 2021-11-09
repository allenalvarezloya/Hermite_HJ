% We simulate the one-dimensional Burgers' equation 
% using 80 cells and m = 3 derivatives
for n = 80
    xl = 0.0;                              % Left endpoint
    xr = 2*pi;                             % Right endpoint
    m = 3;                                 % Number of derivatives used 
    h = (xr - xl)/n;                       % Cell width
    x = xl + h*(0:n)';                     % Primal grid
    xd = xl + h*(1:n)' - 0.5*h;            % Dual grid
    xd_ext = [xd(end); xd; xd(1)];         % Dual grid extended (for BCs)
    u = initialData(x,h,m);                % Set initial data
    ud = zeros(2*m+2,n);                   % Allocate for solution on dual grid
    Hmap = Hermite_map(m,0.0,1.0,0.5,0);   % Hermite interpolation map
    t = 0.0;
    T = 1.5;
    while t < T
        dt = 0.05*h;
        if (t + dt > T)
            dt = T - t;
        end
        dth = 0.5*dt;
        Sd = Sensor(2*m+3,m,u,x);
        Vsd = viscosity(Sd);
        lambda = max(abs((1/h)*u(2,:)));
        v0 = 0.5*(lambda*h)/(2*m+1);
        for i = 1:n
            u_int = Hmap*[u(1:m+1,i); u(1:m+1,i+1)];
            ud(:,i) = RK4visc(u_int,v0*Vsd(i),dth,h,m);
        end
        ud_ext = [ud(:,end) ud ud(:,1)]; % Extend grid for periodic BCs
        Sp = Sensor(2*m+3,m,ud_ext,xd_ext);
        Vsp = viscosity(Sp);
        for i = 1:n+1
            u_int = Hmap*[ud_ext(1:m+1,i); ud_ext(1:m+1,i+1)];
            u(:,i) = RK4visc(u_int,v0*Vsp(i),dth,h,m);
        end
        t = t + dt;
        TITLE = sprintf('Burgers Equation time = %3.2f',t);
        plot(x,u(1,:),'LineWidth',2);
        title(TITLE)
        axis([0 2*pi -1.1 1.1])
        legend('Numerical Solution')
        drawnow
    end
end

