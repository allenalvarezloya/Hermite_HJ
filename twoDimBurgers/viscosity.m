function [Vs] = viscosity(s,CASE)
v = zeros(size(s));
Vs = zeros(size(s));
for k = 1:length(s)
    for j = 1:length(s)
        if s(k,j) < 2
            v(k,j) = 0;
        elseif (2 <= s(k,j) && s(k,j) <= 3)
            v(k,j) = 0.5*(1 + sin(pi*(s(k,j) - 2)/(2)));
        elseif 3 < s(k,j)
            v(k,j) = 1;
        end
    end
end
V = 1 - v;
%% Smooth the viscosity %%
% for k = 2:length(s) - 1
%     Vs(k) = 0.25*(V(k-1) + 2*V(k) + V(k+1));
% end
% Vs(1) = 0.25*(V(length(s)) + 2*V(1) + V(2));
% Vs(length(s)) = 0.25*(V(length(s)-1) +2*V(length(s)) + V(1));
for k = 1:1
    Vs = smooth_visc(V,s,CASE);
    V = Vs;
end
Vs = V;
end