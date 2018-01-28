file = fopen('build/output.bin', 'r');
input = fread(file, [2 9], 'double');
fclose(file);

error = abs(input(2, :) - input(2, end));
t = 10.^(0:-1:-8);

loglog(t, error, '--s')
grid on