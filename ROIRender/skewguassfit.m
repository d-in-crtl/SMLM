function [calcu_ydata] = skewguassfit(x0,xdata)

% expression: y0 + a*gauss*cfgauss
% x0 = [a, xc, sigma, alpha, y0]

y0 = x0(5);
alpha = x0(4);
a = x0(1);
sigma = x0(3);
xc = x0(2);

t = (xdata - xc)./sigma;

fgauss = 1./(sqrt(2*pi).*sigma).*exp(-0.5.*t.^2);
cfgauss = (1+erf(alpha./sqrt(2).*t))./2;
calcu_ydata = y0 + a.*fgauss.*cfgauss;