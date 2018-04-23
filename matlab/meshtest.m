clear
clc
close all

fileID = fopen('../output/bnd.bin', 'wb');
N = 20;
M = 40;

[c, c1, c2, c3, c4] = ring(N, M, 1, 2);
% [c, c1, c2, c3, c4] = nozzel(N, M);
plotarrow(gca, c1(1,:), c1(2,:), 1);
plotarrow(gca, c2(1,:), c2(2,:), 2);
plotarrow(gca, c3(1,:), c3(2,:), 3);
plotarrow(gca, c4(1,:), c4(2,:), 4);
axis equal

fwrite(fileID, [0, N, M] , 'int');
fwrite(fileID, c, 'double');
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
