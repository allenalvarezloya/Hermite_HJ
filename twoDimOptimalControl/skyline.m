function [q_bar] = skyline(Coeff,q)
q_bar = zeros(q+1,1);
for n = 1:q
    q_bar(n) = Coeff(n);
    for j = n+1:q+1
        if (abs(Coeff(j)) > abs(q_bar(n)))
            q_bar(n) = Coeff(j);
        end
    end
end
n = q+1;
q_bar(n) = Coeff(n);
j = q;
if (abs(Coeff(j)) > abs(q_bar(n)))
    q_bar(n) = Coeff(j);
end
end