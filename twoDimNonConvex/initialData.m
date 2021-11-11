function u =  initialData(x,y,mx,my,hx,hy)
	u = zeros(mx+1,my+1,length(x),length(y)); 
    for i = 1:length(x)
        for j = 1:length(y)
            omega_x = 1;
            for idx = 1:mx
                scale_x = omega_x^idx*hx^idx/factorial(idx);
                case_x = mod(idx,4);
                if case_x == 0
                    u(1+idx,1,i,j) = scale_x*sin(x(i));
                elseif case_x == 1
                    u(1+idx,1,i,j) = scale_x*cos(x(i));
                elseif case_x == 2
                    u(1+idx,1,i,j) = -scale_x*sin(x(i));
                elseif case_x == 3
                    u(1+idx,1,i,j) = -scale_x*cos(x(i));
                end
            end
            omega_y = 1;
            for idy = 1:my
                scale_y = omega_y^idy*hy^idy/factorial(idy);
                case_y = mod(idy,4);
                if case_y == 0
                    u(1,1+idy,i,j) = scale_y*cos(y(j));
                elseif case_y == 1
                    u(1,1+idy,i,j) = -scale_y*sin(y(j));
                elseif case_y == 2
                    u(1,1+idy,i,j) = -scale_y*cos(y(j));
                elseif case_y == 3
                    u(1,1+idy,i,j) = scale_y*sin(y(j));
                end
            end

            u(1,1,i,j) = sin(x(i)) + cos(y(j));
        end
    end

end