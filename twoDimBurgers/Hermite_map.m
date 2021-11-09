function tmat =  Hermite_map(m,xl,xr,xc,icase)

% ! This subroutine computes the coefficient matrix tmat which
% ! transfers derivative data of order through m at xl and xr
% ! to the function and derivative data at xc
% ! That is, if p is a polynomial of degree 2m+1
% !
% !  h^k D^k p (xc)/k! = sum_(j=0)^m tmat(k,j) h^j D^j p(xl)/j!
% !                                + tmat(k,j+m+1) h^j D^j p(xr)/j!
% !
% !  icase < 0 => xc=xl (left boundary case)
% !  icase = 0 => xl < xc < xr
% !  icase > 0 => xc=xr (right boundary case)
% IMPLICIT NONE
% INTEGER, INTENT(IN) :: m,icase
% DOUBLE PRECISION, INTENT(IN) :: xl,xr,xc
% DOUBLE PRECISION, DIMENSION(0:2*m+1,0:2*m+1), INTENT(OUT) :: tmat
% DOUBLE PRECISION :: h,z,zc,adl,adr,sign,c1l,c1r,c2l,c2r
% INTEGER :: i,j,k

% ! Compute in normalized coordinates

tmat = zeros(2*m+2,2*m+2);
h = xr-xl;
z = (xc-xl)/h;
zc = z-1.d0;

if (icase > 0)
  z=1.d0;
  zc=0.d0;
elseif (icase < 0)
  z=0.d0;
  zc=-1.d0;
end

bcofs = binomial(m+1);

% !
% ! We begin with the Hermite-Lagrange interpolants:
% !
% !   Q_j (z) = z^(m+1) sum_{k=j}^m h_{kj} (z-1)^k,
% !
% !   j=0, ... , m
% !
% !   satisfying Q_j = (z-1)^j + O((z-1)^(m+1)) at z=1
% !
% !   After some algebra once can show:
% !
% !   h_jj = 1,  h_kj = -sum_{p=j}^{k-1} b_(k-p)^(m+1) h_pj ,
% !              for k>j
% !
% !   here b_(k-p)^(m+1) is the binomial coefficient (m+1)!/((k-p)!(m+1-k+p)!)
% !
% ! To construct the matrix we
% ! now evaluate the interpolant and its derivatives at z
% !
% ! Data to the left is handled by a simple change of variables z -> 1-z
% !
% ! Start with the last column - note that we directly use the recursive
% ! definition of the h's to fold previously computed columns into old
% ! ones. Note that the polynomial is now centered about the midpoint
% !

for i=0:2*m+1
  adl=0.d0;
  adr=0.d0;
  for j=max(0,i-m):min(i,m+1)
    if ((m-i+j) == 0 )
      c2l = 1.d0;
      c2r = 1.d0;
    else
      c2l = z^(m-i+j);
      c2r = zc^(m-i+j);
    end
    if ((m+1-j) == 0)
      c1l = 1.d0;
      c1r = 1.d0;
    else
      c1l = (zc^(m+1-j));
      c1r = (z^(m+1-j));
    end
    adr=adr+bcofs(1+m+1,1+j)*bcofs(1+m,1+i-j)*c1r*c2r;
    adl=adl+bcofs(1+m+1,1+j)*bcofs(1+m,1+i-j)*c1l*c2l;
  end
  tmat(1+i,1+2*m+1)=adr;
  tmat(1+i,1+m)=((-1.d0)^(m+1))*adl;
end

% ! Now loop over the other columns backwards

for k=m-1:-1:0
  for i=0:2*m+1
    adl=0.d0;
    adr=0.d0;
    for j=max(0,i-k):min(i,m+1)
      if ((k-i+j) == 0 )
        c2l = 1.d0;
        c2r = 1.d0;
      else
        c2l = z^(k-i+j);
        c2r = zc^(k-i+j);
      end
      if ((m+1-j) == 0)
        c1l = 1.d0;
        c1r = 1.d0;
      else
        c1l = (zc^(m+1-j));
        c1r = (z^(m+1-j));
      end
      adr=adr+bcofs(1+m+1,1+j)*bcofs(1+k,1+i-j)*c1r*c2r;
      adl=adl+bcofs(1+m+1,1+j)*bcofs(1+k,1+i-j)*c1l*c2l;
    end
    tmat(1+i,1+k+m+1)=adr;
    tmat(1+i,1+k)=((-1.d0)^(m+1))*adl;
    sign=1.d0;
    for j=k+1:m
      sign=-sign;
      tmat(1+i,1+k)=tmat(1+i,1+k)-sign*bcofs(1+m+1,1+j-k)*tmat(1+i,1+j);
      tmat(1+i,1+k+m+1)=tmat(1+i,1+k+m+1)-bcofs(1+m+1,1+j-k)*tmat(1+i,1+j+m+1);
    end
  end
end


function coeffs = binomial(m)
% ! Computes the binomial coefficients of order up through m
% IMPLICIT NONE
% INTEGER, INTENT(IN) :: m
% DOUBLE PRECISION, DIMENSION(0:m,0:m), INTENT(OUT) :: coeffs
% INTEGER :: i,j

coeffs = zeros(m+1,m+1);
coeffs(1,1)=1.d0;

for i=1:m
  coeffs(1+i,1+0)=1.d0;
  for j=1:i
    coeffs(1+i,1+j)=coeffs(1+i,1+j-1)*(i-j+1)/(j);
  end
end
