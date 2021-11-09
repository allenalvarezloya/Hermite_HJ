function [S] = Sensor(q,F,x)
[nodes, weights, P] = lgl_nodes_weights(q);
qh = (q+1)/2;
sol_mat_left = zeros(qh,q+1);
sol_mat_right = zeros(qh,q+1);
scaled_nodes = 0.5*nodes;
scaled_nodes_left = scaled_nodes(qh+1:end)+0.5;
scaled_nodes_right = scaled_nodes(1:qh)-0.5;
for j = 0:q
    sol_mat_left(:,j+1) = scaled_nodes_left.^j;
    sol_mat_right(:,j+1) = scaled_nodes_right.^j;
end
S = zeros(length(x)-1,1);
for k = 1:length(x) - 1
    u_left = sol_mat_left*F(:,k);
    u_right = sol_mat_right*F(:,k+1);
    f = [u_right; u_left];
    C = Coeffs(f,q,P,weights);

    N = length(C)+1;
    b = zeros(length(C),1);
    for n = 1:N-1
        b(n) = 1/n^N;
    end
    b = b/sqrt(sum(b.^2));
    C = sqrt(C.^2 + max(sum(C.^2),1e-5)*b.^2);
    C = skyline(C,q);
    S(k) = Smoothness_param(C,q);

end