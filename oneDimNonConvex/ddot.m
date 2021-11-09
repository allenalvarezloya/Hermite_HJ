function [OUT] = ddot(V,U,k)
  OUT=0;
  for j=1:k-1
    OUT=OUT+j*U(j+1)*V(k-j+1);
  end
  OUT=OUT/k;
end