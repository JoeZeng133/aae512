function [xc, yc, N, M] = readmesh(filename)
%READMESH read mesh from binary file

fileID = fopen(filename, 'rb');
N = fread(fileID, 1, 'int');
M = fread(fileID, 1, 'int');
xc = fread(fileID, M * N, 'double');
yc = fread(fileID, M * N, 'double');
xc = reshape(xc, [M N])';
yc = reshape(yc, [M N])';
fclose(fileID);
end

