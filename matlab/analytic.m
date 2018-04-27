function Fn = analytic(T0, T1, T2, a, rspan, tspan, ax)
% calculate the analytic solution by expanding the solution in Bessel's
% functions

r1 = rspan(1);
r2 = rspan(end);

Tinf = @(r) (T2 - T1) / (log(r2) - log(r1)) * (log(r) - log(r1)) + T1;

eigen_eq = @(l) besselj(0, l * r2) .* bessely(0, l * r1) - besselj(0, l * r1) .* bessely(0, l * r2);
lspan = linspace(0, 500, 6000);
y0 = abs(eigen_eq(lspan));

x0 = diff(y0) > 0;
x1 = diff(y0) < 0;
x1 = [x0, false] & [false, x1];

if ax ~= 0
    plot(lspan, y0, 'linewidth', 1.5), hold on
    plot(lspan(x1), y0(x1), 'sq')
    xlabel('\lambda [m^{-1}]')
    ylabel('|R_0(\lambda)|')
end


l0 = lspan(x1);
l1 = zeros(size(l0));
num_l = numel(l1);
for i = 1 : num_l
    l1(i) = fzero(eigen_eq, l0(i));
end

R0 = @(r, l) besselj(0, l * r) .* bessely(0, l * r1) - besselj(0, l * r1) .* bessely(0, l * r);

Fn = zeros([numel(tspan) numel(rspan)]);

for i = 1 : num_l
    arg1 = integral(@(r) R0(r, l1(i)) .* (T0 - Tinf(r)) .* r, r1, r2);
    Fn = Fn + exp(-a * l1(i)^2 * tspan)' * pi^2 / 2 * l1(i)^2 * besselj(0, l1(i) * r2)^2 / (besselj(0, l1(i) * r1)^2 - besselj(0, l1(i) * r2)^2) * ...
        R0(rspan, l1(i)) * arg1;
end
Fn = Fn + Tinf(rspan);
Fn = Fn';
end


