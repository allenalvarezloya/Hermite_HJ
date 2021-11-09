function [urefined] = refineSol(u,mx,my,nx,ny,refine)
refineGrid = -0.5 + (1/refine)*(0:refine);
u1 = eval_Hermite(u,mx,my,1,1,refineGrid(1),refineGrid(1));
urefined_x = [];
utemp = zeros(refine);
ulast = zeros(refine+1,refine);
urefined_temp = [];
urefined = [];
for j_index = 1:ny
    for i_index = 1:nx
        for i = 1:refine
            for j = 1:refine
                utemp(i,j) = eval_Hermite(u,mx,my,i_index,j_index,refineGrid(i),refineGrid(j));
            end
        end
        urefined_x = [urefined_x; utemp];
    end
    for i = 1:refine+1
        for j = 1:refine
            ulast(i,j) = eval_Hermite(u,mx,my,nx+1,j_index,refineGrid(i),refineGrid(j));
        end
    end
    urefined_temp = [urefined_x; ulast];
    urefined = [urefined urefined_temp];
    urefined_x = [];
    utemp = zeros(refine);
    ulast = zeros(refine+1,refine);
    urefined_temp = [];
end
for i_index = 1:nx
    for i = 1:refine
        for j = 1:refine+1
            utemp(i,j) = eval_Hermite(u,mx,my,i_index,ny+1,refineGrid(i),refineGrid(j));
        end
    end
    urefined_x = [urefined_x; utemp];
end
for i = 1:refine+1
    for j = 1:refine+1
        ulast(i,j) = eval_Hermite(u,mx,my,nx+1,ny+1,refineGrid(i),refineGrid(j));
    end
end
urefined_temp = [urefined_x; ulast];
urefined = [urefined urefined_temp];
urefined(1:refine/2,:) = [];
urefined(:,1:refine/2) = [];
urefined(end-(refine/2-1):end,:) = [];
urefined(:,end-(refine/2-1):end) = [];
