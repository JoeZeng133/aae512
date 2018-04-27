clear
clc
% close all

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
% 
% figure(1)
% for i = 1 : N
%     plot(x(i, :), y(i, :), 'r'), hold on
% end
% 
% for i = 1 : M
%     plot(x(:, i), y(:, i), 'b'), hold on
% end
% axis equal
% xlabel('x [m]')
% ylabel('y [m]')
% fileID = fopen('../output/mesh.bin', 'wb');
% fwrite(fileID, [N, M], 'int');
% fwrite(fileID, x', 'double');
% fwrite(fileID, y', 'double');
% fclose(fileID);

[c, c1, c2, c3, c4] = ring(N, M, 0.1, 0.3, 1, 1.02);
fileID = fopen('../output/bnd.bin', 'wb');
fwrite(fileID, [1, N, M] , 'int');
fwrite(fileID, c, 'double');
fclose(fileID);
!./mesh

fileID = fopen('../output/config.txt', 'wb');
fprintf(fileID, '%.5e ', [T0 T1 T2 a dt]);
fprintf(fileID, '%d', P);
fclose(fileID);

tspan = (1:P) * dt;
rspan = linspace(ri, ro, 1000);
% figure
Fa = analytic(T0, T1, T2, a, rspan, tspan, 0);
%semi-infinite approximation, works for t < 0.3 at least
Fs = T1 + (T0 - T1) * erf((rspan - ri)' * 1 ./ sqrt(4 * a * tspan));
Ftot = [Fs(:, tspan <= 0.3) Fa(:, tspan > 0.3)];

% figure(2)
% plot(rspan, Fa(:, 1),'.-'), hold on
% plot(rspan, Fa(:, 10),'-')
% plot(rspan, Fa(:, 30),'--'), hold off
% xlabel('x [m]')
% ylabel('T [^0C]')
% legend({['t = ', num2str(tspan(1)), 's'], ['t = ', num2str(tspan(10)), 's'], ['t = ', num2str(tspan(30)), 's']}, 'fontsize', 15)

% figure(3)
% slice = find(rspan < 0.14);
% slice = slice(1:5:end);
% 
% plot(rspan(slice), Fa(slice, 20), 'r-'), hold on
% plot(rspan(slice), Fs(slice, 20),'bsq')
% 
% plot(rspan(slice), Fa(slice, 100),'r-')
% plot(rspan(slice), Fs(slice, 100), 'bsq')
% 
% plot(rspan(slice), Fa(slice, 1000), 'r-')
% plot(rspan(slice), Fs(slice, 1000), 'bsq')
% xlabel('x [m]')
% ylabel('T [^0C]')
% legend({'SOP solution', 'Semi-infinite solution'}, 'fontsize', 15)
% 

% !./mesh

%%
figure
[xc, yc] = mesh_plot('../output/mesh.bin', gca);
r = sqrt(xc.^2 + yc.^2);

%%
fileID = fopen('../output/res.bin', 'rb');
data = fread(fileID, N * M * P, 'double');
data = reshape(data, [N M P]);
fclose(fileID);

figure(2)
tp = 10;
pos = rspan < 0.12;
plot(rspan(pos), Ftot(pos, tp), 'linewidth', 1.5), hold on

pos2 = r(N/2,:) < 0.12;
plot(r(N/2,pos2), data(N/2, pos2, tp), 'sq-')

% figure
% plot(rspan, Ftot(:, 100), 'r', 'linewidth', 1.5), hold on
% plot(r(1, :), data(1, :, 100), 'bsq', 'linewidth', 1.5)
% plot(rspan, Ftot(:, 1000), 'r', 'linewidth', 1.5)
% plot(r(1, :), data(1, :, 1000), 'bsq', 'linewidth', 1.5)
% legend({'Analytic solution', 'Numerical solution'}, 'Fontsize', 15)
% xlabel('r [m]')
% ylabel('T [^0C]')
% annotation('textarrow',[0.448214285714286 0.391071428571429],...
%     [0.544238095238096 0.490476190476191],'String',{'t = 100 s'});
% annotation('textarrow',[0.357142857142857 0.298214285714285],...
%     [0.302380952380953 0.294238095238096],'String',{'t = 10 s'});


