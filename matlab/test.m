close all
clc
clear

d = ptsclustering(40, 1, 1.1);
plot(d, zeros(size(d)), 'sq')
min(d)
max(d)