function u =  initialData(x,y,mx,my,hx,hy)
	u = zeros(mx+1,my+1,length(x),length(y)); 
	omega = 1;
    for i = 1:length(x)
        for j = 1:length(y)
            for idy = 0:my
                for idx = 0:mx 
                    scale = omega^(idx+idy)*hx^idx*hy^idy/(factorial(idx)*factorial(idy));
                    CASE = mod(idx + idy,4);
                    if CASE == 0
                        u(1 + idx,1 + idy,i,j) = -scale*cos(x(i)+y(j));
                    elseif CASE == 1
                        u(1 + idx,1 + idy,i,j) = scale*sin(x(i)+y(j));
                    elseif CASE == 2
                        u(1 + idx,1 + idy,i,j) = scale*cos(x(i)+y(j));
                    elseif CASE == 3
                        u(1 + idx,1 + idy,i,j) = -scale*sin(x(i)+y(j));
                    end
                end
            end
        end
    end

end