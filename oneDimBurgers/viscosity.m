function [Vs] = viscosity(s)
v = zeros(length(s),1);
for k = 1:length(s)
    if s(k) < 2
        v(k) = 0;
    elseif (2 <= s(k) && s(k) <= 3)
        v(k) = 0.5*(1 + sin(pi*(s(k) - 2)/(2)));
    elseif 3 < s(k)
        v(k) = 1;
    end
end
V = 1 - v;
%% Smooth the viscosity 
Vs = smooth_visc(V,s);

end