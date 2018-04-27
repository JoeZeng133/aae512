function [xc, yc] = plotmesh(xc, yc, ax)
%MESH_PLOT plot mesh using mesh point coordinates
%   Detailed explanation goes here
N = size(xc, 1);
M = size(xc, 2);
for i = 1 : N
    plot(ax, xc(i, :), yc(i, :), 'r'), hold on
end
for i = 1 : M
    plot(ax, xc(:, i), yc(:, i), 'b'), hold on
end
axis equal
xlabel('x [m]')
ylabel('y [m]')
end

