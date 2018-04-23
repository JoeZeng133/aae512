clear
clc

ri = 0.1;
ro = 0.3;
N = 50;
M = 100;

e1 = linspace(0, 1, N);
e2 = linspace(0, 1, M);

[E1, E2] = ndgrid(e1, e2);
r = E2 * (ro - ri) + ri;
th = E1 * 2 * pi;

x = r .* cos(th);
y = r .* sin(th);

figure(1)
for i = 1 : N
    plot(x(i, :), y(i, :), 'r'), hold on
end

for i = 1 : M
    plot(x(:, i), y(:, i), 'b'), hold on
end

axis equal

fileID = fopen('../output/meshtest.bin', 'wb');
fwrite(fileID, [N, M], 'int');
fwrite(fileID, x', 'double');
fwrite(fileID, y', 'double');
fclose(fileID);


%%
fileID = fopen('../output/res.bin', 'rb');
data = fread(fileID, [N M], 'double');
fclose(fileID);

figure(2)
plot(r(1,:), data(1,:))