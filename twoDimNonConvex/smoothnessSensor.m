function [s] = smoothnessSensor(mx,my,u,nx,ny)
s = zeros(nx,ny);
q = 2*mx+1;
qh = (q+1)/2;
[nodes, weights, P] = lgl_nodes_weights(q);
l_nodes = (nodes(1:qh)+1)/2;
r_nodes = (nodes(qh+1:q+1)-1)/2;
u_bl = zeros(qh,qh);
u_br = zeros(qh,qh);
u_tl = zeros(qh,qh);
u_tr = zeros(qh,qh);
gamma_inv = (2*(0:q)+1)/2;
B = zeros(q+1,q+1);
for y_index = 1:ny
    for x_index = 1:nx
        for j = 1:qh
            for i = 1:qh
                u_bl(i,j) = eval_Hermite(u,mx,my,x_index,y_index,l_nodes(i),l_nodes(j));
                u_br(i,j) = eval_Hermite(u,mx,my,x_index+1,y_index,r_nodes(i),l_nodes(j));
                u_tl(i,j) = eval_Hermite(u,mx,my,x_index,y_index+1,l_nodes(i),r_nodes(j));
                u_tr(i,j) = eval_Hermite(u,mx,my,x_index+1,y_index+1,r_nodes(i),r_nodes(j));
            end
        end
        u_S = [u_bl u_tl; u_br u_tr]; % Check this 
        
        for k = 0:q
            g = sum(weights.*u_S.*P(:,1+k),1);
            for l = 0:q 
                B(1+k,1+l) = sum(weights.*g'.*P(:,1+l))*gamma_inv(1+k)*gamma_inv(1+l);
            end
        end
        b_ordered = zeros(q,1);
        temp = 0.0;
        % This groups coefficients with the same total degree
        for i = 1:q
            for j = 0:i
            b_ordered(i) = max(temp,abs(B(1+i-j,1+j)))+10e-16;
            temp = b_ordered(i);
            end
            temp = 0.0;
        end  
        % Baseline
        Baseline = zeros(length(b_ordered),1);
        N = length(b_ordered)+1;
        for n = 1:N-1
            Baseline(n) = 1/n^N;
        end
        Baseline = Baseline/sqrt(sum(Baseline.^2));
        b_ordered = sqrt(b_ordered.^2 + max(sum(b_ordered.^2),1e-5)*Baseline.^2);
        b_ordered = skyline(b_ordered,q-1);
        s(x_index,y_index) = Smoothness_param(b_ordered,q-1);
    end
end

end