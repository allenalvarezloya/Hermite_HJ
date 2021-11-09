function [PDE] = PDEvisc(u,nu,h,m,Cx,Sx,Sy)
ux = (1/h)*[(1:2*m+1)'.*u(2:2*m+2,:);zeros(1,2*m+2)];
uy = (1/h)*[(1:2*m+1).*u(:,2:2*m+2) zeros(2*m+2,1)];
% uxx = (1/h^2)*[(1:2*m)'.*(2:2*m+1)'.*u(3:2*m+2,:); zeros(2,2*m+2)];
% uyy = (1/h^2)*[(1:2*m).*(2:2*m+1).*u(:,3:2*m+2) zeros(2*m+2,2)];
Sy_uxLong = conv2(Sy,ux);
Sy_ux = Sy_uxLong(1:2*m+2,1:2*m+2);
sign_uy = zeros(2*m+2);
sign_uy(1,1) = sign(uy(1,1));
Sx_sign_uy_uyLong = conv2((Sx+sign_uy),uy);
Sx_sign_uy_uy = Sx_sign_uy_uyLong(1:2*m+2,1:2*m+2);
Sy2Long = conv2(Sy,Sy);
Sy2 = Sy2Long(1:2*m+2,1:2*m+2);
ONE = zeros(2*m+2,2*m+2);
ONE(1,1) = 1;


PDE = ONE - Cx + 0.5*Sy2 - Sx_sign_uy_uy - Sy_ux;
end