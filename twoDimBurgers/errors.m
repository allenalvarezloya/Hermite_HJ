% Run this script to compute the errors in the 
% L1, L2 and Linf norms along with estimated rates 
% of convergence 
Ncells = [20 40 80];
L1error = zeros(length(Ncells),1);
L2error = zeros(length(Ncells),1);
Linferror = zeros(length(Ncells),1);
Count = 1;
for N = Ncells
    m = 3;              % Number of derivatives
    mx = m;             % Number of derivatives in x
    my = m;             % Number of derivatives in y
    xl = 0.0;           % Left endpoint
    xr = 2*pi;          % Right endpoint
    yb = 0.0;           % Lower endpoint
    yt = 2*pi;          % Upper endpoint

    nx = N;             % Number of cells in x
    ny = N;             % Number of cells in y
    hx = (xr-xl)/nx;    % Cell width in x direction
    hy = (yt-yb)/ny;    % Cell width in y direction
    x = xl+(0:nx)'*hx;  % x grid 
    y = yb+(0:ny)'*hy;  % y grid
    Smooth_dual = zeros(N,N);
    Smoooth_primal = zeros(N,N);
    u = initialData(x,y,2*mx+1,2*my+1,hx,hy);
    ud = zeros(2*mx+2,2*my+2,nx,ny);
    ud_extended = zeros(2*mx+2,2*my+2,nx+2,ny+2);
    Hmap = Hermite_map(m,0.0,1.0,0.5,0);
    t = 0.0;
    T = 0.1;
    while t < T
        dt = 0.01*hx;
        if (t + dt > T)
            dt = T - t;
        end
        dth = 0.5*dt;
        Smooth_dual = smoothnessSensor(mx,my,u,nx,ny);
        Vsd = viscosity(Smooth_dual,1);
        lambda_x = max(max(1/hx*abs(u(2,1,:,:))));
        lambda_y = max(max(1/hy*abs(u(1,2,:,:))));
        lambda = max(max(abs(1/hx*u(2,1,:,:)+1/hy*u(1,2,:,:))));
        v0 = 3*lambda*hx/(2*mx+1);
        for i = 1:nx
            for j = 1:ny
                u_int = interp_2D(u(1:mx+1,1:my+1,i,j),u(1:mx+1,1:my+1,1+i,j),...
                    u(1:mx+1,1:my+1,1+i,1+j),u(1:mx+1,1:my+1,i,1+j),m,Hmap);
                ud(:,:,i,j) = RK4_visc(u_int,v0*Vsd(i,j),dth,hx,m);
            end
        end
        % Extend the grid for periodic boundary conditions
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
                u(:,:,1+i,1+j) = RK4_visc(u_int,v0*Vsp(1+i,1+j),dth,hx,m);
            end
        end
        t = t + dt
        u_plot = zeros(nx+1,ny+1);
        for i = 1:nx+1
            for j = 1:ny+1
                u_plot(i,j) = u(1,1,i,j);
            end
        end
%         TITLE = sprintf('T = %3.3f',t);
%         surf(Smooth_primal)
%         title(TITLE)
%         view([40 40])
%         shading interp
%         drawnow
    end
    u_refined = refineSol(u,mx,my,nx,ny,4);
    NN = length(u_refined)-1;
    hh = (2*pi/NN);
    HH = (2*pi/(2*NN));
    XX = 0 + HH*(0:2*NN);
    u_actual = zeros(NN+1,1);
    options = optimset('TolX',1.e-14);
    % Actual solution using Lax-Hopf formula 
    for k = 0:2*NN
        f = @(y) (XX(1+k) - y)^2/(2*T) - cos(2*y);
        [minimizer, fmin] = fminbnd(f,-1000000,1000000,options);
        u_actual(1+k) = fmin;
    end

    U_actual = zeros(NN+1);
    for i = 0:NN
        for j = 0:NN
            U_actual(1+i,1+j) = u_actual(1+i+j);
        end
    end

    L1error(Count) = hh^2*sum(sum(abs(U_actual - u_refined)));
    L2error(Count) = sqrt(sum(sum((U_actual - u_refined).^2)*hh^2));
    Linferror(Count) = max(max(U_actual - u_refined));
    Count = Count +1;
end

conv1 = log2(L1error(1:end-1)./L1error(2:end));
conv2 = log2(L2error(1:end-1)./L2error(2:end));
convinf = log2(Linferror(1:end-1)./Linferror(2:end));

fprintf("\\begin{table}[ht]\n ")
fprintf(sprintf("\\\\caption{Example 1 mx = %i} \n",m))
fprintf("\\begin{center} \n ")
fprintf("\\begin{tabular}{|c|c|c|c|c|c|c|} \n")
fprintf("\\hline \n")
fprintf(" n & $L_1$ error & Convergence & $L_2$ error & Convergence & $L_{\\infty}$ error & Convergence \\\\ \n ")
fprintf("\\hline \n")
fprintf(sprintf(' %i & %3.2e & - & %3.2e & - & %3.2e & - \\\\\\\\ \n',Ncells(1),L1error(1),L2error(1),Linferror(1)))
for z = 2:length(L1error)
    fprintf(sprintf('%i & %3.2e & %3.2f & %3.2e & %3.2f & %3.2e & %3.2f \\\\\\\\ \n',Ncells(z),L1error(z),conv1(z-1),L2error(z),conv2(z-1),Linferror(z),convinf(z-1)))
end
fprintf("\\hline \n")
fprintf("\\end{tabular} \n")
fprintf("\\end{center} \n")
fprintf("\\end{table} \n")
