function [PDE] = PDEvisc(u,nu,h,m)
ux = (1/h)*[(1:2*m+1)'.*u(2:2*m+2,:);zeros(1,2*m+2)];
uy = (1/h)*[(1:2*m+1).*u(:,2:2*m+2) zeros(2*m+2,1)];
uxx = (1/h^2)*[(1:2*m)'.*(2:2*m+1)'.*u(3:2*m+2,:); zeros(2,2*m+2)];
uyy = (1/h^2)*[(1:2*m).*(2:2*m+1).*u(:,3:2*m+2) zeros(2*m+2,2)];
PDE_long = conv2(ux,uy);

PDE = -PDE_long(1:2*m+2,1:2*m+2) + nu*(uxx+uyy);
end