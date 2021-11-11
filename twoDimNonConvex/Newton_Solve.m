function X_sol = Newton_Solve(a,b,x,y,t)
    tol = 10e-16;
    count = 1;
    while sqrt((x - a - t*sin(y))^2 + (y - b + t*cos(x))^2) > tol
        J = [1 -t*cos(y); -t*sin(x) 1];
        F = [x - a - t*sin(y); y - b + t*cos(x)];
        dx = -J\F;
        X = [x;y] + dx;
        x = X(1);
        y = X(2);
        count = count+1;
        if count > 1000
            break
        end
    end
    X_sol = [x;y];
end