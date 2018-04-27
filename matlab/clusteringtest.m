clear
clc
close all

T0 = 0;
T1 = 100;
T2 = 0;
ri = 0.1;
ro = 0.3;
a = 143 / (2.8e3 * 795);

dt = 0.01;
N = 80;
M = 40;
P = 1000;

% [x, y, r] = alge_mesh(N, M, ri, ro, 0.5, 1.02);

% first clustering
c = ring(N, M, 0.1, 0.3, 1, 1.1);
fileID = fopen('../output/bnd1.bin', 'wb');
fwrite(fileID, [1, N, M] , 'int');
fwrite(fileID, c, 'double');
fclose(fileID);
% finer clustering
c = ring(N, M, 0.1, 0.3, 1, 1.01);
fileID = fopen('../output/bnd2.bin', 'wb');
fwrite(fileID, [1, N, M] , 'int');
fwrite(fileID, c, 'double');
fclose(fileID);

!./main ../output/bnd1.bin ../output/mesh1.bin
!./main ../output/bnd2.bin ../output/mesh2.bin
[x1, y1] = readmesh('../output/mesh1.bin');
r1 = sqrt(x1.^2 + y1.^2);
figure(1)
plotmesh(x1, y1, gca);
[x2, y2] = readmesh('../output/mesh2.bin');
r2 = sqrt(x2.^2 + y2.^2);
figure(2)
plotmesh(x2, y2, gca);

%%
% writing configuration
fileID = fopen('../output/config.txt', 'wb');
fprintf(fileID, '%.5e ', [T0 T1 T2 a dt]);
fprintf(fileID, '%d', P);
fclose(fileID);


% solving
!./main ../output/mesh1.bin ../output/config.txt ../output/res1.bin
!./main ../output/mesh2.bin ../output/config.txt ../output/res2.bin
fileID = fopen('../output/res1.bin', 'rb');
data1 = fread(fileID, N * M * P, 'double');
data1 = reshape(data1, [N M P]);
fclose(fileID);

fileID = fopen('../output/res2.bin', 'rb');
data2 = fread(fileID, N * M * P, 'double');
data2 = reshape(data2, [N M P]);
fclose(fileID);

% analytic solution
tspan = (1:P) * dt;
rspan = linspace(ri, ro, 1000);
Fa = analytic(T0, T1, T2, a, rspan, tspan, 0);
Fs = T1 + (T0 - T1) * erf((rspan - ri)' * 1 ./ sqrt(4 * a * tspan));
Ftot = [Fs(:, tspan <= 0.3) Fa(:, tspan > 0.3)];

%% plotting figures
figure(3)
tp = 10;
slice = r1(N/2,:) < 0.12;
plot(r1(N/2,slice), data1(N/2, slice, tp), 'rsq', 'markersize', 10), hold on

slice = r2(N/2,:) < 0.12;
plot(r2(N/2,slice), data2(N/2, slice, tp), 'bo', 'markersize', 10)

slice = rspan < 0.12;
plot(rspan(slice), Ftot(slice, tp), 'k','linewidth', 1.5)

xlabel('x [m]')
ylabel('T [0^C]')
legend({'Coarser clustering', 'Refined clustering', 'Analytic Solution'}, 'fontsize', 15)