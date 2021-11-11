function [u] = interp_2D(BL,BR,TR,TL,m,Hmap)
    u = zeros(2*m+2);
    % Interpolate with respect to x first and then y
    Bottom_int = zeros(2*m+2,m+1);
    Top_int = zeros(2*m+2,m+1);
    for k = 0:m
        Bottom_int(:,1+k) = Hmap*[BL(:,1+k); BR(:,1+k)];
        Top_int(:,1+k) = Hmap*[TL(:,1+k); TR(:,1+k)];
    end
    for k = 0:2*m+1
        u(1+k,:) = Hmap*[Bottom_int(1+k,:)'; Top_int(1+k,:)'];
    end

end