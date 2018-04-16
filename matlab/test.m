clear
clc
close all

fileID = fopen('../output/bnd.bin', 'wb');
N = 20;
M = 20;
fwrite(fileID, N, 'int');
fwrite(fileID, M, 'int');

Mpts = linspace(0, 1, M);
Npts = linspace(0, 1, N);

% exponential clustering
% r = 4;
% Mpos = -1 / (1 - exp(-r)) * (exp(-r * Mpts) - 1);
% Npos = -1 / (1 - exp(-r)) * (exp(-r * Npts) - 1);
% Mpos = fliplr(1 - Mpos);
% Npos = fliplr(1 - Npos);

% power clustering
r = 2;
Mpos = linspace(0, 1, M).^r;
Npos = linspace(0, 1, N).^r;

ang = (sin(linspace(-pi/2, pi/2, N)) + 1) / 2 * pi/2;


c1 = [1 + 0.4 * (1 - fliplr(Mpos));zeros([1 M])];
c2 = [1.4 * cos(ang).^(2/3);sin(ang).^(2/3)];
c3 = [zeros([1 M]); 1 - Mpos];
c4 = [linspace(0, 1, N);zeros([1 N])];

plotarrow(gca, c1(1,:), c1(2,:));
plotarrow(gca, c2(1,:), c2(2,:));
plotarrow(gca, c3(1,:), c3(2,:));
plotarrow(gca, c4(1,:), c4(2,:));
axis equal

fwrite(fileID, c1(:,1:end-1), 'double');
fwrite(fileID, c2(:,1:end-1), 'double');
fwrite(fileID, c3(:,1:end-1), 'double');
fwrite(fileID, c4(:,1:end-1), 'double');
fclose(fileID);

%%

fileID = fopen('../output/mesh.bin', 'rb');
xc = fread(fileID, N * M, 'double');
yc = fread(fileID, N * M, 'double');

xc = reshape(xc, [M N]);
yc = reshape(yc, [M N]);

for i = 1 : M
    plot(xc(i, :), yc(i, :), 'r'), hold on
end

for i = 1 : N
    plot(xc(:, i), yc(:, i), 'r'), hold on
end

axis equal
