function plotarrow(ax, x, y)
%PLOTARROW Summary of this function goes here
%   Detailed explanation goes here
cidx = floor(numel(x) / 2);
cx = x(cidx);
cy = y(cidx);

d = [0;0];
d(1) = x(cidx + 1) - x(cidx - 1);
d(2) = y(cidx + 1) - y(cidx - 1);

d = d / norm(d, 2) * 0.1;

plot(ax, x, y, '*'), hold on
q = quiver(ax, cx, cy, d(1), d(2), 'r');
q.MaxHeadSize = 2;
end

