function [bins, Norm, PairCorrelation] = smPairCorrelation_CPP(coor_1, coor_2, roiszx, roiszy)

% coors are [x, y] in the unit of nm
% roisz is in the unit of nm
% output is 3D with the indexing of [r_12, r_13, delta_theta] (col major)

x_1 = single(coor_1(:, 1));
y_1 = single(coor_1(:, 2));
x_2 = single(coor_2(:, 1));
y_2 = single(coor_2(:, 2));

roisz = roiszx * roiszy;
density_1 = size(coor_1, 1)/roisz;
density_2 = size(coor_2, 1)/roisz;

[bins, Norm, corr] = smPairCorrCPP(x_1, y_1, x_2, y_2, roiszx, roiszy);

Norm = (reshape(Norm, [], size(coor_1, 1)))';
corr = (reshape(corr, [], size(coor_1, 1)))';
PairCorrelation = (sum(corr, 1))';
PairCorrelation = PairCorrelation ./(density_1*density_2) / roisz;




