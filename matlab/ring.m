function [ c, c1, c2, c3, c4 ] = ring(N, M, ri, ro)
%NOZZEL Summary of this function goes here
%   Detailed explanation goes here
% 
% ri = 1;
% ro = 2;

Mpos = linspace(0, 1, M);

% exponential clustering
% r = 4;
% Npts = linspace(0, 1, N);
% Npos = -1 / (1 - exp(-r)) * (exp(-r * Npts) - 1);
% Npos = fliplr(1 - Npos);
% linear
Npos = linspace(0, 1, N);


ang = linspace(0, 2 * pi, M);


c1 = [ri + (ro - ri) * Npos; zeros([1 N])];
c2 = [ro * cos(ang);ro * sin(ang)];
c3 = fliplr(c1);
c4 = fliplr([ri * cos(ang);ri * sin(ang)]);

c = [c1(:,1:end-1), c2(:, 1:end-1), c3(:, 1:end - 1), c4(:, 1:end - 1)];
end

