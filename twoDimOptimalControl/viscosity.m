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
% for k = 1:3
%     Vs = smooth_visc(V,s,CASE);
%     V = Vs;
% end
Vs = V;
end