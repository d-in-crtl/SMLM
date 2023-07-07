function [out] = movie_correction(in, P, Q)
% Transplant of the POLY_2D function from IDL, which works on image while this one only work for given coordinates in 'in'
% INPUT:  in = img to correct
%         P  = array given by IDL polywarp function (the first 16 numbers (1:16) in the .map file given by channel_map.pro)
%         Q  = array given by IDL polywarp function (the second 16 numbers (17:32) in the .map file given by channel_map.pro)
% OUTPUT: out = image corrected

% Algorithm: please refer to the POLY_2D function in IDL
degree = sqrt(size(P, 1)) - 1;

[sizey, sizex] = size(in);
out = zeros(sizey, sizex, 'uint16');

x_in = (1:sizex)' - 1;
y_in = (1:sizey)' - 1;

[XX_in, YY_in] = meshgrid(x_in, y_in);
xx_in = reshape(XX_in, [], 1);
yy_in = reshape(YY_in, [], 1);

xx_out = zeros(size(xx_in ,1), 1);
yy_out = zeros(size(yy_in, 1), 1);

for idx = 1:(degree+1)^2
    pow_x = floor((idx-1)/(degree+1));
    pow_y = mod((idx-1), (degree+1));
    xx_out = xx_out + round(P(idx).*xx_in.^pow_x.*yy_in.^pow_y);
    yy_out = yy_out + round(Q(idx).*xx_in.^pow_x.*yy_in.^pow_y);
end

mask = xx_out < 0 | xx_out > sizex-1 | yy_out < 0 | yy_out > sizey-1;
xx_in(mask) = [];
yy_in(mask) = [];
xx_out(mask) = [];
yy_out(mask) = [];

idx_in = xx_in.*sizey + yy_in + 1;
idx_out = xx_out.*sizey + yy_out + 1;

out(idx_in) = in(idx_out);

