%% demonstrating mesh generation of a nozzel geometry
close all
clc
clear

if ismac
    exec = './main';
elseif ispc
    exec = 'main.exe';
else
    disp('Platform Not Supported')
    return
end

N = 20;
M = 40;
b = 1.01;
[ c, c1, c2, c3, c4 ] = nozzel(N, M, b);
fileID = fopen('../output/bnd.bin', 'wb');
fwrite(fileID, [1, N, M] , 'int');
fwrite(fileID, c, 'double');
fclose(fileID);

system([exec, ' ../output/bnd.bin ../output/mesh.bin']);
[x, y] = readmesh('../output/mesh.bin');
plotmesh(x, y, gca);

