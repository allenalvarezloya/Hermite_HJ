function [s] = sensor(mx,my,u,nx,ny)
q = 2*mx+1;
s = zeros(nx,ny);
[nodes, weights, P] = lgl_nodes_weights(q);
left_nodes = (nodes(1:(q+1)/2)+1)/2
right_nodes = (nodes((q+1)/2+1:q-1)/2

s = 1;
end