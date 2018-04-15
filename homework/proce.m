% Appendix B: Plot Generator
clear
clc
close all

file = fopen('build/output.bin', 'r');
input = fread(file, [2 9], 'double');
fclose(file);

error = abs(input(2, :) - input(2, end));
t = 10.^(0:-1:-8);

loglog(t, error, '--s', 'linewidth', 2)
% axis equal
grid on

xlabel('Step')
ylabel('u')
set(gca, 'FontSize', 15)



polyfit(log10(t(1:end - 1)), log10(error(1: end - 1)), 1)
