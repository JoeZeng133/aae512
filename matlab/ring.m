function [ c, c1, c2, c3, c4 ] = ring(N, M, ri, ro, a, b)
%NOZZEL Summary of this function goes here
%   Detailed explanation goes here
% 

rpos = ptsclustering(M, a, b);
angpos = linspace(0, 1, N);
ang = 2 * pi * angpos;


c1 = fliplr([ri * cos(ang);ri * sin(ang)]);
c2 = [ri + (ro - ri) * rpos; zeros(size(rpos))];
c3 = [ro * cos(ang);ro * sin(ang)];
c4 = fliplr([ri + (ro - ri) * rpos; zeros(size(rpos))]);

c = [c1(:,1:end-1), c2(:, 1:end-1), c3(:, 1:end - 1), c4(:, 1:end - 1)];
end

