clear
clc

t0 = 0;
t1 = 100;
t2 = 20;
r1 = 0.1;
r2 = 0.3;
a = 143 / (2.8e3 * 795);

Tinf = @(r) (t2 - t1) / (log(r2) - log(r1)) * (log(r) - log(r1)) + t1;

eigen_eq = @(l) besselj(0, l * r2) .* bessely(0, l * r1) - besselj(0, l * r1) .* bessely(0, l * r2);
lspan = linspace(0, 500, 6000);
y0 = abs(eigen_eq(lspan));

x0 = diff(y0) > 0;
x1 = diff(y0) < 0;
x1 = [x0, false] & [false, x1];

figure
plot(lspan, y0), hold on
plot(lspan(x1), y0(x1), 'sq')

l0 = lspan(x1);
l1 = zeros(size(l0));
num_l = numel(l1);
for i = 1 : num_l
    l1(i) = fzero(eigen_eq, l0(i));
end


%%
clc

eigen_func = @(r, l) besselj(0, l * r) .* bessely(0, l * r1) - besselj(0, l * r1) .* bessely(0, l * r);
c = zeros(size([1 num_l]));

tstep = 100;
rstep = 1000;

rspan = linspace(r1, r2, rstep);
tspan = linspace(0, 100, tstep);

Fn = zeros([tstep rstep]);

for i = 1 : num_l
    arg1 = integral(@(r) eigen_func(r, l1(i)) .* (t0 - Tinf(r)) .* r, r1, r2);
    Fn = Fn + exp(-a * l1(i)^2 * tspan)' * pi^2 / 2 * l1(i)^2 * besselj(0, l1(i) * r2)^2 / (besselj(0, l1(i) * r1)^2 - besselj(0, l1(i) * r2)^2) * ...
        eigen_func(rspan, l1(i)) * arg1;
end

Fn = Fn + Tinf(rspan);
figure
plot(rspan, Fn(end,:))


%% animation
% interval = 0.01;
% tic
% ymin = min(Fn(:));
% ymax = max(Fn(:));
% 
% for i = 1 : tstep
%     t1 = toc;
%     
%     plot(rspan, Fn(i, :))
%     axis([-inf inf ymin ymax])
%     while(toc - t1 < interval)
%     end
%     getframe;
% end



