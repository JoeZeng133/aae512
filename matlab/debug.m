clear
clc

ri = 1;
ro = 2;

M = 100;
N = 50;

e1 = linspace(0, 1, N);
e2 = linspace(0, 1, M);

[E1, E2] = ndgrid(e1, e2);
r = E2 * (ro - ri) + ri;
th = E1 * 2 * pi;

x = r .* cos(th);
y = r .* sin(th);

% x = E1.^2;
% y = E2.^2;
% 
% J11 = 0.5 ./ sqrt(x);
% J12 = 0;
% J21 = 0;
% J22 = 0.5 ./ sqrt(y);
% L1 = -0.25 ./ x.^1.5;
% L2 = -0.25 ./ y.^1.5;

B = [-1 0.5;-2 6];
A = inv(B);
% 
% x = B(1,1) * E1.^2 + B(1,2) * E2.^2;
% y = B(2,1) * E1.^2 + B(2,2) * E2.^2;


arg1 = A(1,1) * x + A(1,2) * y;
arg2 = A(2,1) * x + A(2,2) * y;
J11 = 0.5 ./ sqrt(arg1) * A(1,1);
J12 = 0.5 ./ sqrt(arg1) * A(1,2);
J21 = 0.5 ./ sqrt(arg2) * A(2,1);
J22 = 0.5 ./ sqrt(arg2) * A(2,2);

J = J11 .* J22 - J21 .* J12;
L1 = -0.25 * (A(1,1)^2 + A(1,2)^2) ./ arg1.^1.5;
L2 = -0.25 * (A(2,1)^2 + A(2,2)^2) ./ arg2.^1.5;

fileID = fopen('../output/meshtest.bin', 'wb');
fwrite(fileID, [N, M], 'int');
fwrite(fileID, x', 'double');
fwrite(fileID, y', 'double');
fclose(fileID);



g11_a = J11.^2 + J12.^2;
g22_a = J21.^2 + J22.^2;
g12_a = J11 .* J21 + J22 .* J12;

g11_a = g11_a(2:end-1, 2:end-1);
g22_a = g22_a(2:end-1, 2:end-1);
g12_a = g12_a(2:end-1, 2:end-1);
L1 = L1(2:end-1, 2:end-1);
L2 = L2(2:end-1, 2:end-1);
J = J(2:end-1, 2:end-1);

%%

fileID = fopen('../output/debug.txt', 'r');

L1_n = fscanf(fileID, '%f', [N-2,M-2]);
L2_n = fscanf(fileID, '%f', [N-2,M-2]);

g11_n = fscanf(fileID, '%f', [N-2,M-2]);
g22_n = fscanf(fileID, '%f', [N-2,M-2]);
g12_n = fscanf(fileID, '%f', [N-2,M-2]);


norm(abs((L1_n - L1) ./ L1), inf)
norm(abs((L2_n - L2) ./ L2), inf)
norm(abs((g11_n - g11_a) ./ g11_a), inf)
norm(abs((g22_n - g22_a) ./ g22_a), inf)
norm(abs((g12_n - g12_a) ./ g12_a), inf)

fclose(fileID);
