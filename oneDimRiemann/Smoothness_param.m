function [s] = Smoothness_param(Coeff,q)
% Compute the matrix for least squares
A = zeros(q,2);
A(:,1) = ones(q,1);
for n = 1:q
    A(n,2) = -log10(n);
end
mapped_q_bar = log10(abs(Coeff)+1e-14);
Y = mapped_q_bar(2:q+1);

% Least Squares
%[param, ~, ~] = housels(A,Y);
param = A\Y;
s = param(2);

end