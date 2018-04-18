function [ c, c1, c2, c3, c4 ] = nozzel(N, M)
%NOZZEL Summary of this function goes here
%   Detailed explanation goes here

Mpos = linspace(0, 1, M);
% exponential clustering
r = 4;
Npts = linspace(0, 1, N);
Npos = -1 / (1 - exp(-r)) * (exp(-r * Npts) - 1);
Npos = fliplr(1 - Npos);


ang = (sin(linspace(-pi/2, pi/2, M)) + 1) / 2 * pi/2;


c1 = [1 + 0.4 * (1 - fliplr(Npos));zeros([1 N])];
c2 = [1.4 * cos(ang).^(2/3);sin(ang).^(2/3)];
c3 = [zeros([1 N]); 1 - Npos];
c4 = [1 * Mpos;zeros([1 M])];

c = [c1(:,1:end-1), c2(:, 1:end-1), c3(:, 1:end - 1), c4(:, 1:end - 1)];
end

