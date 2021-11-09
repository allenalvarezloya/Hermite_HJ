function [u] = initialData(x,h,m)
    % Compute the degrees of freedom for sin(x)
    n = length(x);
    u = zeros(2*m+2,n);
    for i = 1:n
        for idx = 0:2*m+1
            scale = h^idx/factorial(idx);
            CASE = mod(idx,4);
            if(CASE == 0)
                u(1+idx,i) = scale*sin(x(i));
            elseif(CASE == 1)
                u(1+idx,i) = scale*cos(x(i));
            elseif(CASE == 2)
                u(1+idx,i) = -scale*sin(x(i));
            elseif(CASE == 3)
                u(1+idx,i) = -scale*cos(x(i));
            end
        end
    end
    
    
end