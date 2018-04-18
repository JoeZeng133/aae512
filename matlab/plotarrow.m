function plotarrow(ax, x, y, num)
%PLOTARROW Summary of this function goes here
%   Detailed explanation goes here
cx = x(1);
cy = y(1);

d = [0;0];
d(1) = x(2) - x(1);
d(2) = y(2) - y(1);

d = d / norm(d, 2) * 0.1;

plot(ax, x, y, '-+', 'LineWidth', 1.5), hold on
q = quiver(ax, cx, cy, d(1), d(2), 'r');
q.MaxHeadSize = 2;
q.LineWidth = 2;
q.Marker = 'o';

text(ax, x(floor(end/2)), y(floor(end/2)), num2str(num), 'FontSize', 15);
end

