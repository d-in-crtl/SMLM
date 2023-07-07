function K = getBSplineKernel(knlsz, knlorder, wvorder, wvscale)
% get a kernel to smooth image using a B-spline wavelet kernel 
% INPUT:
%   knlsz: scaler, size of the kernel;
%   knlorder: scaler, returns the constructed b-spline wavelet if knlorder == 1, higher order will be inserting zeros in the wavelet according to the a trous algorithm
%   wvorder: the order of the b-spline wavelet, default = 3
%   wvscale: the wavelet will be expanded by wvscale times, default = 2
% OUTPUT:
%   K: a matlab structure that contains upto knlorder constructed wvlet

if nargin < 3
    wvorder = 3;
    wvscale = 2;
elseif nargin < 4
    wvscale = 2;
end

knots = (0 : wvorder);
knot_cntr = wvorder / 2; 

pp = bspline(knots);

sample = (1 : knlsz) / wvscale;
sample_cntr = (sample(1) + sample(end)) / 2;
sample = sample - sample_cntr + knot_cntr;

result = zeros(knlorder, knlsz);

[~, ~ ,bin] = histcounts(sample, knots);
for ii= 1 : length(sample)
    if bin(ii) ~= 0
        pos = bin(ii);
        para = pp.coefs(pos, :);
        for jj = 1 : wvorder
            result(1, ii) = result(1, ii) + para(jj)*(sample(ii) - knots(pos))^(wvorder-jj);
        end
    else
        continue;
    end
end

result(1, :) = result(1, :) / sum(result(1, :));

K = struct([]);
for ii = 1 : knlorder
    K(ii).kernel = reshape(result(1:ii, :), [], 1);
    if ii > 1
        K(ii).kernel(end - (ii-2) : end) = [];
    end
end

