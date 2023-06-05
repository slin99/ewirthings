n = 99;
t_max = 10;
h = 1/(n + 1);
tau = 1/20;

filename = 'WelleNeumann3.gif';
delay = 0.05;

%n=n+2 needed for the generation of T
n = n+2;
j = 1;

T_mesh = 0:tau:t_max;
[X,Y] = meshgrid(0:h:1);

q   = zeros(n^2, t_max/tau + 1);
p   = zeros(n^2, t_max/tau + 1);

f = @(t, x, y) 300 * exp(-10 * t) .* exp(-200 * (x - 0.5).^2 + (y - 0.5).^2);


TT = - sparse(diag(ones(n-1, 1), 1) + diag(ones(n - 1, 1), -1));
TT(1, 2    ) = -2;
TT(n, n - 1) = -2;
T = TT + 4 * speye(n,n);

A = - (kron(speye(n,n), T) + kron(TT, speye(n,n))) / (h * h);
A_plus  = speye(n^2, n^2) + (tau / 2)^2 * A;
A_minus = speye(n^2, n^2) - (tau / 2)^2 * A;


qSurf = surf(X, Y, zeros(n)); zlim([0, 54])
title(gca, 'Wave')
frame = getframe(gcf);
im = frame2im(frame);
[C, map] = rgb2ind(im, 256);
imwrite(C, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', delay);

for t = T_mesh
    r = (f(t + tau, X, Y)' + f(t, X, Y)') / 2;

    q(:, j + 1) = A_minus \ (A_plus * q(:, j) + tau * p(:, j) + tau^2 / 2 * r(:));
    p(:, j + 1) = A_minus \ (A_plus * p(:, j) + tau * A * q(:, j) + tau * r(:));
    
    j = j + 1;
    
    set(qSurf, 'ZData', reshape(q(:, j), [n, n])'); zlim([0, 54])
    title(gca, 'Wave')
    frame = getframe(gcf);
    im = frame2im(frame);
    [C, map] = rgb2ind(im, 256);
    imwrite(C, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', delay);
end
close;