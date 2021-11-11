function [u_val] = eval_Hermite(u,mx,my,x_index,y_index,xloc,yloc)
    u_val = 0.0;
    dfy = 1;
    for j = 1:2*my+2
        dfx = 1;
        for i = 1:2*mx+2
            u_val = u_val+ u(i,j,x_index,y_index)*dfx*dfy;
            dfx = dfx*xloc;
        end
        dfy = dfy*yloc;
    end
end