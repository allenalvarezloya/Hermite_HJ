function [Vs] = viscosity(s)
v = zeros(length(s),1);
V = zeros(length(s),1);
Vs = zeros(length(s),1);
for k = 1:length(s)
    if s(k) < 3
        v(k) = 0;
    elseif (3 <= s(k) && s(k) <= 4)
        v(k) = 0.5*(1 + sin(pi*(s(k) - 3))/(2));
    elseif 4 < s(k)
        v(k) = 1;
    end
end
V = 1 - v;
%% Smooth the viscosity %%
for k = 2:length(s) - 1
    Vs(k) = 0.25*(V(k-1) + 2*V(k) + V(k+1));
end
Vs(1) = 0.25*(V(length(s)) + 2*V(1) + V(2));
Vs(length(s)) = 0.25*(V(length(s)-1) +2*V(length(s)) + V(1));


end