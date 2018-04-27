function writemesh(filename, xc, yc)
%READMESH read mesh from binary file
%   Detailed explanation goes here
fileID = fopen(filename, 'wb');
fwrite(fileID, size(xc), 'int');
fwrite(fileID, xc', 'double');
fwrite(fileID, yc', 'double');
fclose(fileID);
end

