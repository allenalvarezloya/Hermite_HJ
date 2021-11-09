function [u] = initialData(x,h,m)
    n = length(x);
    u = zeros(2*m+2,n);
    for i = 1:n
        omega = pi;
        for idx = 0:2*m+1
            scale = omega^idx*h^idx/factorial(idx);
            CASE = mod(idx,4);
            if(CASE == 0)
                u(1+idx,i) = -scale*cos(omega*x(i));
            elseif(CASE == 1)
                u(1+idx,i) = scale*sin(omega*x(i));
            elseif(CASE == 2)
                u(1+idx,i) = scale*cos(omega*x(i));
            elseif(CASE == 3)
                u(1+idx,i) = -scale*sin(omega*x(i));
            end
        end
    end
    
    
end