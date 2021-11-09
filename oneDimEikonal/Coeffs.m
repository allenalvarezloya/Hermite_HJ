function C = Coeffs(f,q,P,w)
% Analytic form of the inner product squared
gamma_inv = 0.5*(2*(0:q)'+ 1);
% Compute the inner product of the function with 
% Legendre polynomials 
b = zeros(q+1,1);
for k = 1:q+1
    b(k) = sum(f.*w.*P(:,k));
end

% Compute the coefficients by scaling by the norm
C = b.*gamma_inv;
end