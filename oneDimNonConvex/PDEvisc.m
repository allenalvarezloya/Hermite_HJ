function [out] = PDEvisc(u,nu,h,m)
% Compute Derivatives
ux = (1/h)*[(1:2*m+1)'.*u(2:2*m+2); 0];
uxx = (1/h)^2*[(1:2*m)'.*(2:2*m+1)'.*u(3:2*m+2);0;0];
U=[ux(1)+1; ux(2:2*m+2)];
% Apply Taylor recursions
N=2*m+2;
C=zeros(N,1);
S=zeros(N,1);
C(1)=cos(U(1));
S(1)=sin(U(1));
for k=1:N-1
    S(k+1)=C(0+1)*U(k+1)+ddot(C,U,k);
    C(k+1)=-S(0+1)*U(k+1)-ddot(S,U,k);
end
out = C + nu*uxx;

end