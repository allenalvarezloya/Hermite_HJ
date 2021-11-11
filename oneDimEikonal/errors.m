% This script computes the errors in the L-1, L-2 and 
% L-infinity norms along with the estimated rates of convergence
clear, clc;
Ncells = [20 40 80 160];
L1error = zeros(length(Ncells),1);
L2error = zeros(length(Ncells),1);
Linferror = zeros(length(Ncells),1);
count = 1;
for n = Ncells
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
    h_refined = 2*pi/N;
    x_refined = 0 + (0:N)*h_refined;
    u_actual = zeros(N+1,1);
    % Actual solution using Lax-Hopf formula 
    options = optimset('TolX',1.e-14);
    for k = 0:N
        [minimizer, fmin] = fminbnd(@sin,x_refined(k+1)-T,x_refined(k+1)+T,options);
        u_actual(1+k) = fmin;
    end
    L1error(count) = h_refined*sum(abs(u_refined-u_actual));
    L2error(count) = sqrt(h_refined*sum(((u_refined-u_actual).^2)));
    Linferror(count) = max(abs(u_refined - u_actual));
    count = count + 1;
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
for z = 2:length(Ncells)
    fprintf(sprintf('%i & %3.2e & %3.2f & %3.2e & %3.2f & %3.2e & %3.2f \\\\\\\\ \n',Ncells(z),L1error(z),conv1(z-1),L2error(z),conv2(z-1),Linferror(z),convinf(z-1)))
end
fprintf("\\hline \n")
fprintf("\\end{tabular} \n")
fprintf("\\end{center} \n")
fprintf("\\end{table} \n")
