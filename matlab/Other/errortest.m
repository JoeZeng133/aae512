clear
clc
close all

T0 = 0;
T1 = 100;
T2 = 0;
ri = 0.1;
ro = 0.3;
a = 143 / (2.8e3 * 795);

dt = 0.1;
Nlist = [30 40 50 60 70 80];
Mlist = Nlist / 2;
blist = [1.02 1.05 1.08 1.1 1.2 1.5];
num = numel(blist);

data = cell(size(Nlist));
x = cell(size(Nlist));
y = cell(size(Nlist));
r = cell(size(Nlist));
P = 100;

if ismac
    exec = './main';
elseif ispc
    exec = 'main.exe';
else
    disp('Platform Not Supported')
    return
end

% writing configuration
fileID = fopen('../output/config.txt', 'wb');
fprintf(fileID, '%.5e ', [T0 T1 T2 a dt]);
fprintf(fileID, '%d', P);
fclose(fileID);

for i = 1 : num
    N = Nlist(i);
    M = Mlist(i);
    c = ring(N, M, 0.1, 0.3, 1, 5);
    fileID = fopen('../output/bnd.bin', 'wb');
    fwrite(fileID, [1, N, M] , 'int');
    fwrite(fileID, c, 'double');
    fclose(fileID);
    system([exec, ' ../output/bnd.bin ../output/mesh.bin'])
    [x{i}, y{i}] = readmesh('../output/mesh.bin');
    r{i} = sqrt(x{i}.^2 + y{i}.^2);

    system([exec, ' ../output/mesh.bin ../output/config.txt ../output/res.bin'])
    fileID = fopen('../output/res.bin', 'rb');
    data{i} = fread(fileID, N * M * P, 'double');
    data{i} = reshape(data{i}, [N M P]);
    fclose(fileID);
end

% analytic solution
tspan = (1:P) * dt;
rspan = linspace(ri, ro, 1000);
Fa = @(r, t) analytic(T0, T1, T2, a, r, t, 0);
Fs = @(r, t) T1 + (T0 - T1) * erf((r - ri)' * 1 ./ sqrt(4 * a * t));


%%
error = zeros(size(Nlist));
tp = P;
% figure(1)
for i = 1 : num
    xd = r{i}(end/2, :);
    if tp * dt < 0.3
        yn = Fs(xd, tp * dt)';
    else
        yn = Fa(xd, tp * dt)';
    end
    
    yd = data{i}(end/2, :, tp);
    error(i) = max(abs(yd - yn));
end

figure(2)
ele_size = 1 ./ sqrt(Nlist .* Mlist);
p = polyfit(ele_size, error, 2);
plot(ele_size, error, 'sq'), hold on
ele_size_s = linspace(ele_size(1), ele_size(end), 100);
plot(ele_size_s, polyval(p, ele_size_s))
xlabel('Element Size')
ylabel('Error [^0C]')

