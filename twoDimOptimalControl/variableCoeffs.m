function [C,Sx,Sy] =  variableCoeffs(x,y,mx,my,hx,hy)
C = zeros(2*mx+2,2*my+2,length(x),length(y)); 
Sx = zeros(2*mx+2,2*my+2,length(x),length(y)); 
Sy = zeros(2*mx+2,2*my+2,length(x),length(y)); 
omega = 1;
for i = 1:length(x)
    for j = 1:length(y)
        for idx = 0:2*mx+1 
            scale = omega^idx*hx^idx/factorial(idx);
            CASE = mod(idx,4);
            if CASE == 0
                C(1 + idx,1,i,j) = scale*cos(x(i));
            elseif CASE == 1
                C(1 + idx,1,i,j) = -scale*sin(x(i));
            elseif CASE == 2
                C(1 + idx,1,i,j) = -scale*cos(x(i));
            elseif CASE == 3
                C(1 + idx,1,i,j) = scale*sin(x(i));
            end
        end
    end
end

for i = 1:length(x)
    for j = 1:length(y)
        for idx = 0:2*mx+1 
            scale = omega^idx*hx^idx/factorial(idx);
            CASE = mod(idx,4);
            if CASE == 0
                Sx(1 + idx,1,i,j) = scale*sin(x(i));
            elseif CASE == 1
                Sx(1 + idx,1,i,j) = scale*cos(x(i));
            elseif CASE == 2
                Sx(1 + idx,1,i,j) = -scale*sin(x(i));
            elseif CASE == 3
                Sx(1 + idx,1,i,j) = -scale*cos(x(i));
            end
        end
    end
end



for i = 1:length(x)
    for j = 1:length(y)
        for idy = 0:2*my+1 
            scale = omega^idy*hy^idy/factorial(idy);
            CASE = mod(idy,4);
            if CASE == 0
                Sy(1,1 + idy,i,j) = scale*sin(y(j));
            elseif CASE == 1
                Sy(1,1 + idy,i,j) = scale*cos(y(j));
            elseif CASE == 2
                Sy(1,1 + idy,i,j) = -scale*sin(y(j));
            elseif CASE == 3
                Sy(1,1 + idy,i,j) = -scale*cos(y(j));
            end
        end
    end
end


end