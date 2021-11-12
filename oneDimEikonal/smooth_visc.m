function [Vs] = smooth_visc(V,s)
    Vs = zeros(length(s),1);
    for k = 2:length(s) - 1
        Vs(k) = 0.25*(V(k-1) + 2*V(k) + V(k+1));
    end
    Vs(1) = 0.25*(V(length(s)) + 2*V(1) + V(2));
    Vs(length(s)) = 0.25*(V(length(s)-1) +2*V(length(s)) + V(1));
end