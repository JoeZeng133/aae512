function writemesh(filename, xc, yc)
%READMESH write mesh to file

fileID = fopen(filename, 'wb');
fwrite(fileID, size(xc), 'int');
fwrite(fileID, xc', 'double');
fwrite(fileID, yc', 'double');
fclose(fileID);
end

