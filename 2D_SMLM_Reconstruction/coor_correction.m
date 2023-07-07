function [out] = coor_correction(in, P, Q)
% Transplant of the POLY_2D function from IDL, which works on image while this one only work for given coordinates in 'in'
% INPUT:  in = [x, y], NOTE: shouldn't be [y, x]
%         P  = array given by IDL polywarp function (the first 16 numbers (1:16) in the.coor file given by channel_map.pro)
%         Q  = array given by IDL polywarp function (the second 16 numbers (17:32) in the.coor file given by channel_map.pro)
% OUTPUT: out = [x, y], coordinates modified by P and Q

% Algorithm: outx = P(idx).*inx.^(floor(idx/(degree+1))).*iny.^(mod(idx/(degree+1)))

degree = sqrt(size(P, 1)) - 1;

out = in;
in_x = in(:,1);
in_y = in(:,2);

out_x = zeros(size(in_x,1), 1);
out_y = zeros(size(in_y,1), 1);

for idx = 1 : (degree+1)^2
    pow_x = floor((idx-1)/(degree+1));
    pow_y = mod((idx-1), (degree+1));
    out_x = out_x + P(idx).*in_x.^pow_x.*in_y.^pow_y;
    out_y = out_y + Q(idx).*in_x.^pow_x.*in_y.^pow_y;
end
out(:,1) = out_x;
out(:,2) = out_y;
