function [d] = ptsclustering(N, a, b)
%PTSCLUSTERING clusters points at two sides of (0,1)
%  a determines clustering distribution between two sides, b determines
%  level of clustering
x = linspace(0, 1, N);

if a > 0.5
    a1 = 1 - a;
else
    a1 = a;
end

arg1 = (2 * a1 + b) * ((b + 1) / (b - 1)).^((x - a1) / (1 - a1)) + 2 * a1 - b;
arg2 = ((2 * a1 + 1) * (1 + ((b + 1) / (b - 1)).^((x - a1) / (1 - a1))));
d = arg1 ./ arg2;

d = (d - min(d)) / (max(d) - min(d));

if a > 0.5
    d = 1 - fliplr(d);
end

end

