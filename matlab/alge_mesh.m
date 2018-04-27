function [x, y, r] = alge_mesh(N, M, ri, ro, a, b)
%ALGE_MESH Summary of this function goes here
%   Detailed explanation goes here

e1 = linspace(0, 1, N);
e2 = ptsclustering(M, a, b);

[E1, E2] = ndgrid(e1, e2);
r = E2 * (ro - ri) + ri;
th = E1 * 2 * pi;

x = r .* cos(th);
y = r .* sin(th);
end

