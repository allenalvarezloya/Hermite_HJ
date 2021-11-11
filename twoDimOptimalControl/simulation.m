% This script computes and plots the numerical soltuion to 
% a problem related to optimal cost
clear, clc;
N = 40;                                 % Number of Cells
m = 3;                                  % Number of derivatives
mx = m;                                 % Number of derivatives in x
my = m;                                 % Number of derivatives in y
xl = -pi;                               % Left endpoint
xr = pi;                                % Right endpoint
yb = -pi;                               % Lower endpoint
yt = pi;                                % Upper endpoint

nx = N;                                 % Number of cells in x
ny = N;                                 % Number of cells in y
hx = (xr-xl)/nx;                        % Cell width in x direction
hy = (yt-yb)/ny;                        % Cell width in y direction
x = xl+(0:nx)'*hx;                      % x grid 
y = yb+(0:ny)'*hy;                      % y grid
xd = xl+(1:nx)'*hx-0.5*hx;              % x dual grid 
yd = yb+(1:ny)'*hy-0.5*hy;              % y dual grid
Smooth_dual = zeros(N,N);
Smoooth_primal = zeros(N,N);
u = initialData(x,y,mx,my,hx,hy);
ud = zeros(2*mx+2,2*my+2,nx,ny);
ud_extended = zeros(2*mx+2,2*my+2,nx+2,ny+2);
Hmap = Hermite_map(m,0.0,1.0,0.5,0);
% Compute sin and cos Taylor series
Spx = 0*u;
Sdx = 0*ud;
Cpx = 0*u;
Cdx = 0*ud;
Spy = 0*u;
Sdy = 0*ud;
[Cpx,Spx,Spy] = variableCoeffs(x,y,mx,my,hx,hy);
[Cdx,Sdx,Sdy] = variableCoeffs(xd,yd,mx,my,hx,hy);
t = 0.0;
T = 1.0;

while t < T
    dt = 0.1*hx;
    if (t + dt > T)
        dt = T - t;
    end
    dth = 0.5*dt;
    Smooth_dual = smoothnessSensor(mx,my,u,nx,ny);
    Vsd = viscosity(Smooth_dual,1);
    lambda_x = max(sin(y));
    lambda_y = max(cos(x));
    lambda = max(lambda_x,lambda_y);
    v0 = lambda*hx/(2*mx+1);
    for i = 1:nx
        for j = 1:ny
            u_int = interp_2D(u(1:mx+1,1:my+1,i,j),u(1:mx+1,1:my+1,1+i,j),...
                u(1:mx+1,1:my+1,1+i,1+j),u(1:mx+1,1:my+1,i,1+j),m,Hmap);
            ud(:,:,i,j) = RK4_visc(u_int,v0*Vsd(i,j),dth,hx,m,Cdx(:,:,i,j),Sdx(:,:,i,j),Sdy(:,:,i,j));
        end
    end
    ud_extended = 0*ud_extended;
    ud_extended(:,:,2:nx+1,2:ny+1) = ud;
    ud_extended(:,:,1,2:ny+1) = ud(:,:,end,:);
    ud_extended(:,:,end,2:ny+1) = ud(:,:,1,:);
    ud_extended(:,:,2:nx+1,1) = ud(:,:,:,end);
    ud_extended(:,:,2:nx+1,end) = ud(:,:,:,1);
    ud_extended(:,:,1,1) = ud(:,:,end,end);
    ud_extended(:,:,end,1) = ud(:,:,1,end);
    ud_extended(:,:,end,end) = ud(:,:,1,1);
    ud_extended(:,:,1,end) = ud(:,:,end,1);
    Smooth_primal = smoothnessSensor(mx,my,ud_extended,nx+1,ny+1);
    Vsp = viscosity(Smooth_primal,2);
    for i = 0:nx
        for j = 0:ny 
            u_int = interp_2D(ud_extended(1:mx+1,1:my+1,1+i,1+j),ud_extended(1:mx+1,1:my+1,2+i,1+j),...
            ud_extended(1:mx+1,1:my+1,2+i,1+1+j),ud_extended(1:mx+1,1:my+1,1+i,2+j),m,Hmap);
            u(:,:,1+i,1+j) = RK4_visc(u_int,v0*Vsp(1+i,1+j),dth,hx,m,Cpx(:,:,1+i,1+j),Spx(:,:,1+i,1+j),Spy(:,:,1+i,1+j));
        end
    end
    
    t = t + dt;
    u_plot = zeros(nx+1,ny+1);
    for i = 1:nx+1
        for j = 1:ny+1
            u_plot(i,j) = u(1,1,i,j);
        end
    end
    [X,Y] = meshgrid(x,y);
    figure(5)
    TITLE = sprintf('T = %3.3f',t);
    surf(X,Y,u_plot')
    title(TITLE)
    shading interp
    axis([-pi pi -pi pi -1 2.5])
    drawnow
end
u_plot = zeros(nx+1,ny+1);
for i = 1:nx+1
    for j = 1:ny+1
        u_plot(i,j) = sign(u(1,2,i,j));
    end
end
[X,Y] = meshgrid(x,y);
figure(5)
TITLE = sprintf('T = %3.3f',t);
surf(X,Y,u_plot')
shading interp
axis([-pi pi -pi pi -1 2.5])