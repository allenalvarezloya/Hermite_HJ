% We compute the L1, L2 and L infinity errors
% along with the estimated rates of convergence 
% when we have an odd number of cells
Ncells = [41 81 161];
L1error = zeros(length(Ncells),1);
L2error = zeros(length(Ncells),1);
Linferror = zeros(length(Ncells),1);
count = 1;
for n = Ncells
    xl = -1.0;
    xr = 1.0;
    m = 3;
    h = (xr - xl)/n;
    x = xl + h*(0:n)';
    xd = xl + h*(1:n)' - 0.5*h;
    xd_ext = [xd(end); xd; xd(1)];
    ud = zeros(2*m+2,n);
    u = initialData(x,h,m,1);
%     u2 = initialData(x,h,m,2);
%     u = 0.0*u1;
    ux = zeros(2*m+2,n+1);
    udx = zeros(2*m+2,n+2);
    Hmap = Hermite_map(m,0.0,1.0,0.5,0);
    t = 0.0;
    T = 1.0;
    while t < T
        dt = 0.05*h;
        if (t + dt > T)
            dt = T - t;
        end
        dth = 0.5*dt;
        Sd = Sensor(2*m+1,u,x);
        Vsd = viscosity(Sd);
        lambda = 2.0;
        v0 = (lambda*h)/(2*m+1);
        for i = 1:n
            u_int = Hmap*[u(1:m+1,i); u(1:m+1,i+1)];
            ud(:,i) = RK4visc(u_int,v0*Vsd(i),dth,h,m);
        end
%         for i = 0.5*n+1:n
%             u_int = Hmap*[u2(1:m+1,i); u2(1:m+1,i+1)];
%             ud(:,i) = RK4visc(u_int,v0*Vsd(i),dth,h,m);
%         end
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
%         u1 = u;
%         u2 = u;
        t = t + dt;
        TITLE = sprintf('time = %3.2f',t);
        plot(x,u(1,:),'LineWidth',2);
        title(TITLE)
        axis([-1.0 1.0 -2.0 0.0])
        drawnow
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

    options = optimset('TolX',1.e-16);
    rng default % For reproducibility
    opts = optimoptions(@fmincon,'Algorithm','sqp');
    for k = 0:N
        fun = @(v) x_refined(1+k)*v - 0.25*(v.^2-1).*(v.^2-4);
        problem = createOptimProblem('fmincon','objective',...
    fun,'x0',0,'lb',-2,'ub',2,'options',opts);
        gs = GlobalSearch;
        [z,f] = run(gs,problem);
        VALUE = min(f,u_actual(1+k));
        u_actual(1+k) = VALUE;
    end
    L1error(count) = h_refined*sum(abs(u_refined-u_actual));
    L2error(count) = sqrt(h_refined*sum(((u_refined-u_actual).^2)));
    Linferror(count) = max(abs(u_refined - u_actual));
    count = count + 1;
    plot(x_refined(1:10:end),u_refined(1:10:end),'o',x_refined,u_actual,'LineWidth',2)
    axis([-1.0 1.0 -2.0 -1.0])

end
close all
plot(x_refined,u_refined,'LineWidth',2)
axis([-1 1 -2 -0.9])
hold on
plot(x_refined(1:21:end),u_actual(1:21:end),'ro')
hold off
legend('Numerical Solution','Exact Solution')
saveas(gcf,'Example_04_1D','epsc')

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