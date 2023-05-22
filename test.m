function [u_h,err] = leeee(n, w, plot_result)
    h = 1/(n+1);
    omega = 2 * pi * w;

    %n=n+2 needed for the generation of T 
    n = n+2;
    
    [X,Y] = meshgrid(0:h:1);
    
    u       = @(x,y)     sin(omega * (x .* x + y)) + y;
    grad_u  = @(x, y)    [2 * omega *  x .* cos(omega * (x .* x + y)); omega * cos(omega * (x .* x + y))] + [0;1];
    g       = @(n, x, y) dot(grad_u(x, y), n);
    f       = @(x, y)    omega^2 * (4 * x .* x + 1) .* sin(omega * (x .* x + y)) - 2 * omega * cos(omega * (x .* x + y));

    TT = - sparse(diag(ones(n-1, 1), 1) + diag(ones(n - 1, 1), -1));
    TT(1, 2    ) = -2;
    TT(n, n - 1) = -2;
    T = TT + 4 * speye(n,n);

    fg    = f(X,Y);

    left  = arrayfun(@(y) g([1; 0], 0, y), Y(:, 1));
    right = arrayfun(@(y) g([ -1; 0], 1, y), Y(:, 1));
    top   = arrayfun(@(x) g([ 0; 1], x, 1), X(1, :));
    bot   = arrayfun(@(x) g([ 0;-1], x, 0), X(1, :));

    fg(:, 1) = fg(:, 1) + 2 * left  / h;
    fg(:, n) = fg(:, n) + 2 * right / h;
    fg(1, :) = fg(1, :) + 2 * top   / h;
    fg(n, :) = fg(n, :) + 2 * bot   / h;

    fg = fg';
    fh = fg(:);

    %create d for projection
    dg = ones(n,n);
    dg(:,1) = dg(:,1)*1/2;
    dg(:,end) = dg(:,end)*1/2;
    dg(1,:) = dg(1,:)*1/2;
    dg(end,:) = dg(end,:)*1/2;
    d = dg(:);
    
    % todo nicer component wise mult with ones(1,n)*1/2 thus corners -> 1/4
    %Ah = (kron(speye(n,n), T)) + (kron(TT, speye(n,n))) / (h * h);
    Ah = gallery("neumann",n^2)/h^2;
    %disp(sprintf("Mat diff: %f" , norm(Ah-bb,"inf")));
    Ah = [Ah; ones(1, n*n)];

    c = norm(d, 2)^2;
    orthProjector =  @(x) x - dot(d, x) * d / c;
    orthProject=  [orthProjector(fh);0];

    u_h = Ah \ orthProject;
   
    uGrid= u(X,Y);
    uGrid = uGrid';
    uGridVector = uGrid(:);
    errV = orthProjector(uGridVector)-u_h;
    err = norm(errV,"inf");
    
    if plot_result == "true"
        hold on
        tiledlayout(1,2)
        nexttile
        surf(X, Y, reshape(errV,[n,n]));
        title("error")
        nexttile
        surf(X, Y, reshape(u_h, [n, n]));
        title("Solution")
    end
end
