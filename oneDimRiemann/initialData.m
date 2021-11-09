function [u] = initialData(x,h,m,type)
n = length(x);
u = zeros(2*m+2,n);
if type == 1
    for i = 1:n
        if (x(i) > 0)
            u(1+0,i) = -2*x(i);
        u(1+1,i) = -2*h;
        else
            u(1+0,i) = 2*x(i);
        u(1+1,i) = 2*h;
        end
    end
elseif type == 2
    for i = 1:n
        if (x(i) >= 0)
            u(1+0,i) = -2*x(i);
            u(1+1,i) = -2*h;
        else
            u(1+0,i) = 2*x(i);
            u(1+1,i) = 2*h;
        end
    end
end
end