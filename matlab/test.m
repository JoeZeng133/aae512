clear
clc
close all

fileID = fopen('../output/bnd.bin', 'wb');
N = 5;
M = 5;
fwrite(fileID, N, 'int');
fwrite(fileID, M, 'int');

ri = 1;
ro = 2;

ang = linspace(0, 2 * pi, M);
c1 = [linspace(ri, ro, N);zeros([1 N])];
c2 = ro * [cos(ang); sin(ang)];
c3 = fliplr(c1);
c4 = ri * [cos(fliplr(ang)); sin(fliplr(ang))];

plot(c1(1,:),c1(2,:),'-'), hold on
plot(c2(1,:),c2(2,:),'sq')
plot(c3(1,:),c3(2,:),'*')
plot(c4(1,:),c4(2,:))
axis equal

fwrite(fileID, c1(:,1:end-1), 'double');
fwrite(fileID, c2(:,1:end-1), 'double');
fwrite(fileID, c3(:,1:end-1), 'double');
fwrite(fileID, c4(:,1:end-1), 'double');
fclose(fileID);

%%
ksid = 1 / (N - 1);
etad = 1 / (M - 1);

xc = ones([M N]);
yc = ones([M N]);
xc(1, :) = c1(1,:);
xc(:, end) = c2(1,:)';
xc(end,:) = fliplr(c3(1,:));
xc(:, 1) = fliplr(c4(1,:))';

yc(1, :) = c1(2,:);
yc(:, end) = c2(2,:)';
yc(end,:) = fliplr(c3(2,:));
yc(:, 1) = fliplr(c4(2,:))';



[x1, x2] = gradient(xc, ksid, etad);
[y1, y2] = gradient(yc, ksid, etad);

x11 = diff(xc, 2, 2) / ksid^2;
x22 = diff(xc, 2, 1) / etad^2;
y11 = diff(yc, 2, 2) / ksid^2;
y22 = diff(yc, 2, 1) / etad^2;

alpha = x2.^2 + y2.^2;
beta = x1 .* x2 + y1 .* y2;
gamma = x1.^2 + y1.^2;

