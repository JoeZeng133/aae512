function [ c, c1, c2, c3, c4 ] = nozzel(N, M, b)
%NOZZEL produces boundary mesh points for a nozzel geometry
%   Detailed explanation goes here

Mpos = linspace(0, 1, M);
% exponential clustering
Npos = ptsclustering(N, 0, b);


ang = (sin(linspace(-pi/2, pi/2, M)) + 1) / 2 * pi/2;


c1 = [1 + 0.4 * Npos;zeros([1 N])];
c2 = [1.4 * cos(ang).^(2/3);sin(ang).^(2/3)];
c3 = [zeros([1 N]); fliplr(Npos)];
c4 = [Mpos;zeros([1 M])];

c = [c1(:,1:end-1), c2(:, 1:end-1), c3(:, 1:end - 1), c4(:, 1:end - 1)];
end

