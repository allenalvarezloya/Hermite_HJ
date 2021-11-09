function [s] = smoothnessSensor(mx,my,u,nx,ny)
s = zeros(nx,ny);
q = 2*mx+1;
qh = (q+1)/2;
[nodes, weights, P] = lgl_nodes_weights(q);
scaled_nodes = 0.5*nodes;
pos_nodes = scaled_nodes(qh+1:end)+0.5;
neg_nodes = scaled_nodes(1:qh)-0.5;

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
                u_bl(i,j) = eval_Hermite(u,mx,my,x_index,y_index,pos_nodes(i),pos_nodes(j));
                u_br(i,j) = eval_Hermite(u,mx,my,x_index+1,y_index,neg_nodes(i),pos_nodes(j));
                u_tl(i,j) = eval_Hermite(u,mx,my,x_index,y_index+1,pos_nodes(i),neg_nodes(j));
                u_tr(i,j) = eval_Hermite(u,mx,my,x_index+1,y_index+1,neg_nodes(i),neg_nodes(j));
            end
        end
        u_S = [u_tr u_br ; u_tl u_bl];  
        
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
            b_ordered(i) = max(temp,abs(B(1+i-j,1+j)));
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