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
Nlist = [16 20 40 60 80];
Mlist = [8 10 20 30 40];
num = numel(Nlist);
data = cell(size(Nlist));
x = cell(size(Nlist));
y = cell(size(Nlist));
r = cell(size(Nlist));
P = 1000;

% writing configuration
fileID = fopen('../output/config.txt', 'wb');
fprintf(fileID, '%.5e ', [T0 T1 T2 a dt]);
fprintf(fileID, '%d', P);
fclose(fileID);

for i = 1 : numel(Nlist)
    N = Nlist(i);
    M = Mlist(i);
    c = ring(N, M, 0.1, 0.3, 1, 2);
    fileID = fopen('../output/bnd.bin', 'wb');
    fwrite(fileID, [1, N, M] , 'int');
    fwrite(fileID, c, 'double');
    fclose(fileID);
    !./main ../output/bnd.bin ../output/mesh.bin
    [x{i}, y{i}] = readmesh('../output/mesh.bin');
    r{i} = sqrt(x{i}.^2 + y{i}.^2);

    !./main ../output/mesh.bin ../output/config.txt ../output/res.bin
    fileID = fopen('../output/res.bin', 'rb');
    data{i} = fread(fileID, N * M * P, 'double');
    data{i} = reshape(data{i}, [N M P]);
    fclose(fileID);
end

% analytic solution
tspan = (1:P) * dt;
rspan = linspace(ri, ro, 1000);
Fa = analytic(T0, T1, T2, a, rspan, tspan, 0);
Fs = T1 + (T0 - T1) * erf((rspan - ri)' * 1 ./ sqrt(4 * a * tspan));
Ftot = [Fs(:, tspan <= 0.3) Fa(:, tspan > 0.3)];


%%
for i = 1 : num
    xd = r{i}(end/2, :);
    yd = data{i}(end/2, :, 500)' - analytic(T0, T1, T2, a, xd, 500 * dt, 0);
    plot(xd, yd, 'sq-'), hold on
end
