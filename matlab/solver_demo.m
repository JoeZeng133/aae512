%% demonstrate numerical solution and anlytic solution

clear
clc
close all

% determining operating system
if ismac || isunix
    exec = './main';
elseif ispc
    exec = 'main.exe';
else
    disp('Platform Not Supported')
    return
end

% parameters for the system
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

% generating analytic solutions
tspan = (1:P) * dt;
rspan = linspace(ri, ro, 1000);
Fa = analytic(T0, T1, T2, a, rspan, tspan, 0);
Fs = T1 + (T0 - T1) * erf((rspan - ri)' * 1 ./ sqrt(4 * a * tspan));
Ftot = [Fs(:, tspan <= 0.3) Fa(:, tspan > 0.3)]; % for t < 0.1 use semi-infinite solution

% figure that shows oscillations of the SOP solution
figure(1)
plot(rspan, Fa(:, 1),'.-'), hold on
plot(rspan, Fa(:, 10),'-')
plot(rspan, Fa(:, 30),'--'), hold off
xlabel('x [m]')
ylabel('T [^0C]')
legend({['t = ', num2str(tspan(1)), 's'], ['t = ', num2str(tspan(10)), 's'], ['t = ', num2str(tspan(30)), 's']}, 'fontsize', 15)

% figure that shows the agreement between semi-infinite solution and sop
% solution at appropriate time steps
figure(2)
slice = find(rspan < 0.14);
slice = slice(1:5:end);

plot(rspan(slice), Fa(slice, 20), 'r-'), hold on
plot(rspan(slice), Fs(slice, 20),'bsq')

plot(rspan(slice), Fa(slice, 100),'r-')
plot(rspan(slice), Fs(slice, 100), 'bsq')

plot(rspan(slice), Fa(slice, 1000), 'r-')
plot(rspan(slice), Fs(slice, 1000), 'bsq'), hold off
xlabel('x [m]')
ylabel('T [^0C]')
legend({'SOP solution', 'Semi-infinite solution'}, 'fontsize', 15)

%% generating numerical results and compare to analytic one
% generating boundary mesh points for mesh generation
[c, c1, c2, c3, c4] = ring(N, M, 0.1, 0.3, 1, 1.02);
fileID = fopen('../output/bnd.bin', 'wb');
fwrite(fileID, [1, N, M] , 'int');
fwrite(fileID, c, 'double');
fclose(fileID);
system([exec, ' ../output/bnd.bin ../output/mesh.bin']) % run mesh generation program
[xc, yc] = readmesh('../output/mesh.bin');
figure(3)
plotmesh(xc, yc, gca); %mesh plot
r = sqrt(xc.^2 + yc.^2);

% write configuration file for the system
fileID = fopen('../output/config.txt', 'wb');
fprintf(fileID, '%.5e ', [T0 T1 T2 a dt]);
fprintf(fileID, '%d', P);
fclose(fileID);
system([exec, ' ../output/mesh.bin ../output/config.txt ../output/res.bin']); %run solver program

% read result
fileID = fopen('../output/res.bin', 'rb');
data = fread(fileID, N * M * P, 'double');
data = reshape(data, [N M P]);
fclose(fileID);

% figure that shows the analytic and numerical solution
figure(4)
plot(rspan, Ftot(:, 100), 'r', 'linewidth', 1.5), hold on
plot(r(1, :), data(1, :, 100), 'bsq', 'linewidth', 1.5)
plot(rspan, Ftot(:, 1000), 'r', 'linewidth', 1.5)
plot(r(1, :), data(1, :, 1000), 'bsq', 'linewidth', 1.5)
legend({'Analytic solution', 'Numerical solution'}, 'Fontsize', 15)
xlabel('r [m]')
ylabel('T [^0C]')
annotation('textarrow',[0.341071428571427 0.289285714285714],...
    [0.17857142857143 0.202380952380952],'String',{'t = 10 s'});
annotation('textarrow',[0.371428571428571 0.314285714285714],...
    [0.479952380952383 0.426190476190478],'String',{'t = 100 s'});


