function [Correlation, rho_center, rad_center] = smTripleCorrelation_Frc(coor_1, coor_2, coor_3, roiszx, roiszy, MaxRhoSize, rho_res, rad_res)
% coors are [x, y] in the unit of nm
% roisz is in the unit of nm
% MaxRadSize is predefined as 90 while rad_res is predefined as 2*pi/90
% -- NOTE that uint8 can handle MaxRadSize
% ================== Configure MaxRhos and MaxRads ========================

MaxRadSize = ceil(2*pi/rad_res);
rad_res = 2*pi/MaxRadSize;

rho_center = ((1:MaxRhoSize)' - 0.5).*rho_res; 
rad_center = ((1:MaxRadSize)' - 0.5).*rad_res - pi;

% ================== Construct Normalization Base =========================
BaseX_A = bsxfun(@times, rho_center, cos(rad_center'));
BaseY_A = bsxfun(@times, rho_center, sin(rad_center'));

TmpX_A = repmat(BaseX_A', MaxRadSize+1, 1);
TmpY_A = repmat(BaseY_A', MaxRadSize+1, 1);

BaseX_B = zeros(MaxRadSize+1, MaxRadSize*MaxRhoSize);
BaseX_B(:) = TmpX_A; 
BaseX_B(end, :) = [];

BaseY_B = zeros(MaxRadSize+1, MaxRadSize*MaxRhoSize); 
BaseY_B(:) = TmpY_A; 
BaseY_B(end, :) = [];

clear TmpX_A TmpY_A;

% =================== Initialize Correlation Results ======================
Correlation_inner = zeros(MaxRhoSize, MaxRadSize, MaxRhoSize);
Correlation_outer = zeros(MaxRhoSize, MaxRadSize, MaxRhoSize);

NormUnit = bsxfun(@times, rho_center.*rad_res.*rho_res, rho_center'.*rad_res.*rho_res);
NormUnit = repmat(NormUnit, [1, 1, MaxRadSize]);
NormUnit = permute(NormUnit, [1, 3, 2]);

% ====================== Calculating spots density ========================
roisz = roiszx * roiszy;
N1 = size(coor_1, 1);
density_2 = size(coor_2, 1)/roisz;
density_3 = size(coor_3, 1)/roisz;

% ===================== Spliting coordinates 1 ============================
ind_inner = coor_1(:, 1) >= MaxRhoSize*rho_res & coor_1(:, 1) <= roiszx - MaxRhoSize*rho_res & coor_1(:, 2) >= MaxRhoSize*rho_res & coor_1(:, 2) <= roiszy - MaxRhoSize*rho_res;
coor_1_inner = coor_1(ind_inner, :);
coor_1_outer = coor_1(~ind_inner, :);

% ===================== Calculating delta theta and r =====================
deltaX_12_inner = bsxfun(@minus, coor_2(:, 1), coor_1_inner(:, 1)');
deltaY_12_inner = bsxfun(@minus, coor_2(:, 2), coor_1_inner(:, 2)');
[Theta_12_inner, Rho_12_inner] = cart2pol(deltaX_12_inner, deltaY_12_inner);
clear deltaX_12_inner deltaY_12_inner;

deltaX_12_outer = bsxfun(@minus, coor_2(:, 1), coor_1_outer(:, 1)');
deltaY_12_outer = bsxfun(@minus, coor_2(:, 2), coor_1_outer(:, 2)');
[Theta_12_outer, Rho_12_outer] = cart2pol(deltaX_12_outer, deltaY_12_outer);
clear deltaX_12_outer deltaY_12_outer;

deltaX_13_inner = bsxfun(@minus, coor_3(:, 1), coor_1_inner(:, 1)');
deltaY_13_inner = bsxfun(@minus, coor_3(:, 2), coor_1_inner(:, 2)');
[Theta_13_inner, Rho_13_inner] = cart2pol(deltaX_13_inner, deltaY_13_inner);
clear deltaX_13_inner deltaY_13_inner;

deltaX_13_outer = bsxfun(@minus, coor_3(:, 1), coor_1_outer(:, 1)');
deltaY_13_outer = bsxfun(@minus, coor_3(:, 2), coor_1_outer(:, 2)');
[Theta_13_outer, Rho_13_outer] = cart2pol(deltaX_13_outer, deltaY_13_outer);
clear deltaX_13_outer deltaY_13_outer;

% ============= Histogram and Normalize at each coor_1_inner ==============
for ii = 1 : size(coor_1_inner, 1)
    
    Norm = MaxRadSize.*ones(MaxRhoSize, MaxRadSize, MaxRhoSize);
    tmp_Corr = zeros(MaxRhoSize, MaxRadSize, MaxRhoSize);
    
    r_12 = Rho_12_inner(:, ii);
    phi_12 = Theta_12_inner(:, ii);
    ind_outlier_12 = r_12 > MaxRhoSize * rho_res;
    r_12(ind_outlier_12) = [];
    phi_12(ind_outlier_12) = [];
    
    r_13 = Rho_13_inner(:, ii);
    phi_13 = Theta_13_inner(:, ii);
    ind_outlier_13 = r_13 > MaxRhoSize * rho_res;
    r_13(ind_outlier_13) = [];
    phi_13(ind_outlier_13) = [];
    
    delta_PHI = bsxfun(@minus, phi_13, phi_12');
    delta_PHI(delta_PHI > pi) = delta_PHI(delta_PHI > pi) - 2*pi;
    delta_PHI(delta_PHI < -pi) = delta_PHI(delta_PHI < -pi) + 2*pi;
    
    r_13_idx = ceil(r_13 ./ rho_res);
    r_12_idx = ceil(r_12 ./ rho_res);
    dphi_idx = ceil(delta_PHI ./ rad_res) + MaxRadSize / 2;
    
    for mm = 1 : size(r_12_idx, 1)
        for kk = 1 : size(r_13_idx, 1)
            tmp_Corr(r_13_idx(kk), dphi_idx(kk, mm), r_12_idx(mm)) = tmp_Corr(r_13_idx(kk), dphi_idx(kk, mm), r_12_idx(mm)) + 1; % sample size is much smaller than the bin size, so better not use built-in hist functions
        end
    end
    
    tmp_Corr = tmp_Corr ./ (Norm.*NormUnit) ./(density_2*density_3);
    thresh_3_tmp = tmp_Corr(end, :, :);
    thresh_2_tmp = tmp_Corr(:, :, end);
    thresh_3 = mean(thresh_3_tmp(:)) + 3 * std(thresh_3_tmp(:)); 
    thresh_2 = mean(thresh_2_tmp(:)) + 3 * std(thresh_2_tmp(:));
    tmp_ind = tmp_Corr >= thresh_3 & tmp_Corr >= thresh_2;
    tmp_Corr(tmp_ind) = 1;
    tmp_Corr(~tmp_ind) = 0;
    
    Correlation_inner = Correlation_inner + tmp_Corr./N1;

end

% ============= Histogram and Normalize at each coor_1_outer ==============
for ii = 1 : size(coor_1_outer, 1)
    
    % ------------------------ Normalization ------------------------------
    Norm_A = zeros(MaxRhoSize, MaxRadSize); 
    Norm_B = zeros(MaxRadSize, MaxRadSize*MaxRhoSize);
    
    NormX_A = BaseX_A + coor_1_outer(ii, 1);
    NormY_A = BaseY_A + coor_1_outer(ii, 2);
    
    ind_A = NormX_A > 0 & NormX_A < roiszx & NormY_A > 0 & NormY_A < roiszy;
    Norm_A(ind_A) = 1;
    
    NormX_B = BaseX_B + coor_1_outer(ii, 1);
    NormY_B = BaseY_B + coor_1_outer(ii, 2);
    ind_B = NormX_B > 0 & NormX_B < roiszx & NormY_B > 0 & NormY_B < roiszy;
    Norm_B(ind_B) = 1;
    
    Norm = Norm_A * Norm_B;
    Norm = reshape(Norm, MaxRhoSize, MaxRadSize, MaxRhoSize);
    Norm = fftshift(Norm, 2);
    
    ind_zero_mask = Norm == 0;
    
    % ---------------- 3D histogram delta theta and r ---------------------
    tmp_Corr = zeros(MaxRhoSize, MaxRadSize, MaxRhoSize);
    
    r_12 = Rho_12_outer(:, ii);
    phi_12 = Theta_12_outer(:, ii);
    ind_outlier_12 = r_12 > MaxRhoSize * rho_res;
    r_12(ind_outlier_12) = [];
    phi_12(ind_outlier_12) = [];
    
    r_13 = Rho_13_outer(:, ii);
    phi_13 = Theta_13_outer(:, ii);
    ind_outlier_13 = r_13 > MaxRhoSize * rho_res;
    r_13(ind_outlier_13) = [];
    phi_13(ind_outlier_13) = [];
    
    delta_PHI = bsxfun(@minus, phi_13, phi_12');
    delta_PHI(delta_PHI > pi) = delta_PHI(delta_PHI > pi) - 2*pi;
    delta_PHI(delta_PHI < -pi) = delta_PHI(delta_PHI < -pi) + 2*pi;
    
    r_13_idx = ceil(r_13 ./ rho_res);
    r_12_idx = ceil(r_12 ./ rho_res);
    dphi_idx = ceil(delta_PHI ./ rad_res) + MaxRadSize / 2;
    
    for mm = 1 : size(r_12_idx, 1)
        for kk = 1 : size(r_13_idx, 1)
            tmp_Corr(r_13_idx(kk), dphi_idx(kk, mm), r_12_idx(mm)) = tmp_Corr(r_13_idx(kk), dphi_idx(kk, mm), r_12_idx(mm)) + 1;
        end
    end
    
    tmp_Corr = tmp_Corr ./ (Norm.*NormUnit) ./(density_2*density_3);
    tmp_Corr(ind_zero_mask) = 1;
    thresh_3_tmp = tmp_Corr(end, :, :);
    thresh_2_tmp = tmp_Corr(:, :, end);
    thresh_3 = mean(thresh_3_tmp(:)) + 2.25 * std(thresh_3_tmp(:)); 
    thresh_2 = mean(thresh_2_tmp(:)) + 2.25 * std(thresh_2_tmp(:));
    tmp_ind = tmp_Corr >= thresh_3 & tmp_Corr >= thresh_2;
    tmp_Corr(tmp_ind) = 1;
    tmp_Corr(~tmp_ind) = 0;
    
    Correlation_outer = Correlation_outer + tmp_Corr./N1;

end

Correlation = Correlation_outer + Correlation_inner;