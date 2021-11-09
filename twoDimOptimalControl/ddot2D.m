function [h] = ddot2D(v,u,k)
    S = 0.0;
    if k(1) == 0
        kp = k(2);
        for z = 1:kp-1
            S = S + (kp-z)*dot(v(1,1+z),u(1,1+kp-z)); % Fine
        end
    elseif k(2) == 0
        kp = k(1);
        for z = 1:kp-1
            S = S + (kp-z)*dot(v(1+z,1),u(1+kp-z,1)); % Fine
        end
    elseif k(1) == 1 && k(2) == 1
        kp = 1;
        S = S + dot(v(1,2),u(2,1)); % Fine
    else
        kp = min(k(1),k(2));
        if kp == 1 && kp == k(1)
            S = S + dot(v(1,2:k(2)+1),u(2,k(2):-1:1)); % Fine
        elseif kp == 1 && kp == k(2)
            S = S + dot(v(2:k(1)+1,1),u(k(1):-1:1,2)); % Fine
        elseif kp == k(1)
            S = S + kp*dot(v(1,2:k(2)+1),u(1+kp,k(2):-1:1)); % Fine
            for z = 1:kp-1
                S = S + (kp-z)*dot(v(z+1,1:k(2)+1),u(1+kp-z,1+k(2):-1:1)); % Fine
            end
        elseif kp == k(2)
            S = S + kp*dot(v(2:k(1)+1,1),u(k(1):-1:1,1+kp)); % Fine
            for z = 1:kp-1
                S = S + (kp-z)*dot(v(1:k(1)+1,z+1),u(1+k(1):-1:1,1+kp-z)); % Fine
            end
        end
    end
    S = S/kp;
    h = S; 
end
