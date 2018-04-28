clear
clc
close all

T0 = 0;
T1 = 100;
T2 = 20;
ri = 0.1;
ro = 0.3;
a = 143 / (2.8e3 * 795);

dt = 0.1;
N = 50;
M = 25;
P = 1000;

if ismac
    exec = './main';
elseif ispc
    exec = 'main.exe';
else
    disp('Platform Not Supported')
    return
end

% [x, y, r] = alge_mesh(N, M, ri, ro, 0.5, 1.02);
c = ring(N, M, 0.1, 0.3, 1, 2);
fileID = fopen('../output/bnd1.bin', 'wb');
fwrite(fileID, [1, N, M] , 'int');
fwrite(fileID, c, 'double');
fclose(fileID);

system([exec, ' ../output/bnd1.bin ../output/mesh1.bin'])
[x1, y1] = readmesh('../output/mesh1.bin');
r1 = sqrt(x1.^2 + y1.^2);
figure(1)
plotmesh(x1, y1, gca);
%%
% writing configuration
fileID = fopen('../output/config.txt', 'wb');
fprintf(fileID, '%.5e ', [T0 T1 T2 a dt]);
fprintf(fileID, '%d', P);
fclose(fileID);

% solving

system([exec, ' ../output/mesh1.bin ../output/config.txt ../output/res1.bin'])
fileID = fopen('../output/res1.bin', 'rb');
data1 = fread(fileID, N * M * P, 'double');
data1 = reshape(data1, [N M P]);
fclose(fileID);

% analytic solution
tspan = (1:P) * dt;
rspan = linspace(ri, ro, 1000);
Fa = analytic(T0, T1, T2, a, rspan, tspan, 0);
Fs = T1 + (T0 - T1) * erf((rspan - ri)' * 1 ./ sqrt(4 * a * tspan));
Ftot = [Fs(:, tspan <= 0.3) Fa(:, tspan > 0.3)];

%% plotting figures
figure(2)
for tp = [10, 100, 1000]
    plot(rspan, Ftot(:, tp), 'r'), hold on
    plot(r1(N/2,:), data1(N/2,:,tp), 'bsq')
end
xlabel('x [m]')
ylabel('T [0^C]')
legend({'Analytic Solution', 'Numerical Solution'}, 'fontsize', 15)
% Create textarrow
annotation('textarrow',[0.292857142857143 0.330357142857143],...
    [0.164285714285715 0.228571428571429],'String',{'t=1 s'});

% Create textarrow
annotation('textarrow',[0.3625 0.326785714285714],...
    [0.41804761904762 0.36904761904762],'String',{'t = 10 s'});

% Create textarrow
annotation('textarrow',[0.5125 0.4625],...
    [0.597619047619048 0.542857142857143],'String',{'t = 100 s'});


