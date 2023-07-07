function [reports, ims_res, endfrm] = GaussMFA2D_SCMOS_wrapper(ims, startfrm, offset, var, gain, PSFsigma, iterations, boxsz, lim, filtmethod, thresh, algorithm, option)
% ========================================================================= 
% Main function wrapping local maximum (CPU) and sub-pixel localization (GPU)
% Use startfrm and endfrm to prevent from running out of RAM
%
%   INPUTS: 
%       ims:            uint16, 3D raw image stacks
%       startfrm:       index of the first frame w.r.t the raw image stack before calling this function.
%       var:            double, 2D variance map of the used camera chip in ADU^2 (must be the same xy size as im)
%       offset:         double, 2D offsets map of the used camera chip in ADU (must be the same xy size as im)
%       gain:           double, 2D gain map of the used camera chip in ADU/e- (must be the same xy size as im)
%       PSFsigma:       initial guess or pre-calibrated PSFSigma    
%       iterations:     number of iterations for MLE fitting
%       boxsz:          sub_region size
%       lim:            maximum number of subregion concurrently sent to GPU
%       filtmethod:     'box':      use a box kernel for image filtering;
%                       'gauss':    use a gaussian kernel for image filtering;
%                       'wvlet':    use a b-spline wavelet kernel for image filtering;    
%       thresh:         scaler:     constant thresholding (least favoured);
%                       2D matrix:  thresholding using camera pixel-independent std of camera read-out noise
%                       'std':      thresholding using the std of the (raw image - fisrt order smoothing)
%       algorithm:      'NR':       Newton-Raphson method for MLE fitting
%                       'LM':       Levernberg-Marquardt method for MLE fitting
%                       'ABFGSB':   Adapted BFGS method for MLE fitting
%       option:         1:          PSFsigma is fixed during fitting
%                       2:          PSFsigma is not fixed during fitting
%                       3:          PSFsigma along x and y are independent during fitting
%   OUTPUTS: 
%       reports:        fitting result table: [x, y, z, I, bg, num, frm, uncert_x, uncert_y, LLR]
%       endfrm:         index of the last frame w.r.t the raw image stack after calling this function.
% 


%% ============================ Initialization ============================ 
BUFFER = 8E9; % preserve 8G bytes buffer allocated for subregion (float) and subvar (float) 

NMAX = 5; % Up to 5 emitters will be attepmeted for MFA fitting
P_VAL = 0.05; % Creteria determining if n+1 emitters fits better than n emitters 
switch option
    case 1
        VNUM0 = 3; % Number of parameters for each Gaussian emitter
    case 2
        VNUM0 = 4;
    case 3
        VNUM0 = 5;
    otherwise
        error('Please specify the gauss model option (1: fixed PSFsigma, 2: unfixed PSFsigma, 3: independent PSFsigma_x and _y)');
end
LLR_thresh = zeros(NMAX, 1, 'single'); % LLR threshold for push-and-pull during MFA fitting
for i = 1 : NMAX
    LLR_thresh(i) = chi2inv(1 - P_VAL, boxsz * boxsz - VNUM0 * i - 1);
end


% Algorithm Info
switch algorithm
    case 'NR'
        PSF_Fitter = @GPUGaussMFA2D_NR_SCMOS;
    case 'LM'
        PSF_Fitter = @GPUGaussMFA2D_LM_SCMOS;
    case 'ABFGSB'
        PSF_Fitter = @GPUGaussMFA2D_ABFGSB_SCMOS;
    otherwise
        error('Please specify the fitting algorithm ("NR", "LM", "ABFGSB").');
end


%% ============================ Segmentation ==============================
[locations, subregion, subvar, ims_res, endfrm] = ImSeg_SCMOS(BUFFER, ims, startfrm, offset, var, gain, PSFsigma, boxsz, filtmethod, thresh);


%% =========================== MFA MLE Fitting ============================
numOfPeaks = size(locations, 1);

y_top = repelem(locations(:, 1), NMAX);
x_left = repelem(locations(:, 2), NMAX);
frm = repelem(locations(:, 3), NMAX);

tot_sub_x = zeros(NMAX * numOfPeaks, 1);
tot_sub_y = zeros(NMAX * numOfPeaks, 1);
tot_sx = zeros(NMAX * numOfPeaks, 1);
tot_sy = zeros(NMAX * numOfPeaks, 1);
tot_photon = zeros(NMAX * numOfPeaks, 1);

tot_CRLB = zeros(NMAX * numOfPeaks, VNUM0);

tot_bg = zeros(NMAX * numOfPeaks, 1);
tot_num = zeros(NMAX * numOfPeaks, 1);
tot_LLR = zeros(NMAX * numOfPeaks, 1);

segs = ceil(numOfPeaks / lim);
for i = 0 : segs - 1
    
    st_ = i * lim;
    ed_ = min(st_ + lim, numOfPeaks);
    
    switch option
        case 1
            [sub_y, sub_x, sub_photon, sub_bg, CRLB, LLR, nnum] = PSF_Fitter(subregion(:,:,st_+1:ed_), subvar(:,:,st_+1:ed_), PSFsigma, iterations, option, LLR_thresh);
        case 2
            [sub_y, sub_x, sub_sx, sub_photon, sub_bg, CRLB, LLR, nnum] = PSF_Fitter(subregion(:,:,st_+1:ed_), subvar(:,:,st_+1:ed_), PSFsigma, iterations, option, LLR_thresh);
            tot_sx(st_ * NMAX + 1 : ed_ * NMAX) = double(sub_sx);
        case 3
            [sub_y, sub_x, sub_sx, sub_sy, sub_photon, sub_bg, CRLB, LLR, nnum] = PSF_Fitter(subregion(:,:,st_+1:ed_), subvar(:,:,st_+1:ed_), PSFsigma, iterations, option, LLR_thresh);
            tot_sx(st_ * NMAX + 1 : ed_ * NMAX) = double(sub_sx);
            tot_sy(st_ * NMAX + 1 : ed_ * NMAX) = double(sub_sy);
    end

    tot_sub_x(st_ * NMAX + 1 : ed_ * NMAX) = double(sub_x);
    tot_sub_y(st_ * NMAX + 1 : ed_ * NMAX) = double(sub_y);
    tot_photon(st_ * NMAX + 1 : ed_ * NMAX) = double(sub_photon);
    
    tot_CRLB(st_ * NMAX + 1 : ed_ * NMAX, :) = double(CRLB);
    
    tot_bg(st_ * NMAX + 1 : ed_ * NMAX) = repelem(double(sub_bg), NMAX);
    tot_num(st_ * NMAX + 1 : ed_ * NMAX) = repelem(double(nnum), NMAX);
    tot_LLR(st_ * NMAX + 1 : ed_ * NMAX) = repelem(double(LLR), NMAX);
    
    
end

%% ======================= Result Filtering ===============================
mask1 = tot_num == 0 | tot_photon <= eps('single');
mask2 = tot_sub_x < 0 | tot_sub_x > boxsz | tot_sub_y < 0 | tot_sub_y > boxsz;
mask3 = false(numOfPeaks * NMAX, 1);
for i = 1 : VNUM0
    mask3 = mask3 | tot_CRLB(:, i) <= eps('single') | isnan(tot_CRLB(:, i));
end
mask = mask1 | mask2 | mask3;

tot_sub_x(mask) = [];
tot_sub_y(mask) = [];
tot_sx(mask) = [];
tot_sy(mask) = [];
tot_photon(mask) = [];
tot_bg(mask) = [];
tot_num(mask) = [];
tot_LLR(mask) = [];
tot_CRLB(mask, :) = [];

y_top(mask) = [];
x_left(mask) = [];
frm(mask) = [];

reports = [double(x_left) + tot_sub_x, double(y_top) + tot_sub_y, tot_sx, tot_sy, tot_photon, tot_bg, tot_num, double(frm), tot_CRLB(:, 2), tot_CRLB(:, 1), tot_LLR];
