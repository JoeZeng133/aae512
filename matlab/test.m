clear
clc
close all

fileID = fopen('../output/bnd.bin', 'wb');
N = 20;
M = 20;
fwrite(fileID, N, 'int');
fwrite(fileID, M, 'int');

r = 1.5;
Npos = linspace(0, 1, N).^r;
Mpos = linspace(0, 1, M).^r;
ang = linspace(0, pi / 2, N);
c1 = [linspace(0, 0.7, N);zeros([1 N])];
c2 = [0.7 + 0.3 * (1 - fliplr(Mpos));zeros([1 M])];
c3 = [cos(ang);sin(ang)];
c4 = [zeros([1 M]); 1 - Mpos];

plot(c1(1,:),c1(2,:),'*'), hold on
plot(c2(1,:),c2(2,:),'*')
plot(c3(1,:),c3(2,:),'*')
plot(c4(1,:),c4(2,:),'*')
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
