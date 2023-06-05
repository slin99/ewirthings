n = 149;
t_max = 0.25;
h = 1 / (n + 1);
tau = 1e-3;

w = 2.5;
omega = 2 * pi * w;

filename = 'WelleNeumann2.gif';
delay = 0.05;

%n=n+2 needed for the generation of T
n = n+2;
j = 1;

T_mesh = 0:tau:t_max;
[X,Y] = meshgrid(0:h:1);

q   = zeros(n^2, t_max/tau + 1);
p   = zeros(n^2, t_max/tau + 1);
err = zeros(  1, t_max/tau + 1);

u       = @(t, x, y)     cos(4 * pi * t) .* ( sin(omega * (x .* x + y)) + y );
dt_u    = @(t, x, y)     - 4 * pi * sin(4 * pi * t) .* ( sin(omega * (x .* x + y)) + y );
grad_u  = @(t, x, y)     cos(4 * pi * t) .*  [2 * omega * x .* cos(omega * (x .* x + y)); omega * cos(omega * (x .* x + y)) + 1];
g       = @(n, t, x, y)  dot(grad_u(t, x, y), n);
f       = @(t, x, y)     - (4 * pi)^2 * u(t, x, y) + cos(4 * pi * t) .* (omega^2 * (4 * x .* x + 1) .* sin(omega * (x .* x + y)) - 2 * omega * cos(omega * (x .* x + y)));

u0 = u(0, X, Y)';
u1 = dt_u(0, X, Y)';

q(:, 1) = u0(:);
p(:, 1) = u1(:);




TT = - sparse(diag(ones(n-1, 1), 1) + diag(ones(n - 1, 1), -1));
TT(1, 2    ) = -2;
TT(n, n - 1) = -2;
T = TT + 4 * speye(n,n);

A       = - (kron(speye(n), T) + kron(TT, speye(n))) / (h * h);
A_minus = speye(n^2) - tau^2 * A; 
A_plus  = speye(n^2) + tau^2 * A;


qSurf = surf( X, Y, u0'); zlim([-3.5, 3])
title(gca, 'Wave')
frame = getframe(gcf);
im = frame2im(frame);
[C, map] = rgb2ind(im, 256);
imwrite(C, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', delay);




for t = T_mesh
     fg    = f(t, X, Y)';

    left  = arrayfun(@(y) g([-1; 0], t, 0, y), Y(:, 1));
    right = arrayfun(@(y) g([ 1; 0], t, 1, y), Y(:, 1));
    top   = arrayfun(@(x) g([ 0; 1], t, x, 1), X(1, :));
    bot   = arrayfun(@(x) g([ 0;-1], t, x, 0), X(1, :));
    
    fg(:, 1) = fg(:, 1) + 2 * left  / h;
    fg(:, n) = fg(:, n) + 2 * right / h;
    fg(n, :) = fg(n, :) + 2 * top   / h;
    fg(1, :) = fg(1, :) + 2 * bot   / h;

    fg = fg';
    r = fg(:);

    q(:, j + 1) = q(:, j) + tau *     p(:, j);
    p(:, j + 1) = p(:, j) + tau * A * q(:, j) + tau   * r;
    
    err(j) = max(abs(reshape(q(:, j), [n, n])' - u(t, X, Y)), [], 'all');

    j = j + 1;
    
    set(qSurf, 'ZData', reshape(q(:, j), [n, n])'); zlim([-3.5, 3])
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [C, map] = rgb2ind(im, 256);
    imwrite(C, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay);
end

errFig = uifigure('Name', 'WelleNeumann1', 'NumberTitle', 'off');
errAx  = uiaxes(errFig);

plot(errAx, T_mesh, err);
title(errAx, 'Error')

close;