function [locations, subregion, subvar, ims_res, endfrm] = ImSeg_SCMOS(buffer, ims, startfrm, offset, var, gain, PSFsigma, boxsz, filtmethod, thresh)
% ========================================================================= 
% Locate the local maximum of a image stack frame-by-frame
% Crop a square boxsz x boxsz around each local maximum
%
%   INPUTS: 
%       buffer:         byte, number of bytes are reserved for subregions and subvars (RAM protection)
%       ims:            uint16, 3D raw image stacks
%       startfrm:       index of the first frame w.r.t the raw image stack before calling this function.
%       var:            double, 2D variance map of the used camera chip in ADU^2 (must be the same xy size as im)
%       offset:         double, 2D offsets map of the used camera chip in ADU (must be the same xy size as im)
%       gain:           double, 2D gain map of the used camera chip in ADU/e- (must be the same xy size as im)
%       PSFsigma:       initial guess or pre-calibrated PSFSigma    
%       boxsz:          sub_region size
%       filtmethod:     'box':      use a box kernel for image filtering;
%                       'gauss':    use a gaussian kernel for image filtering;
%                       'wvlet':    use a b-spline wavelet kernel for image filtering;    
%       thresh:         scaler:     constant thresholding (least favoured);
%                       2D matrix:  thresholding using camera pixel-independent std of camera read-out noise
%                       'std':      thresholding using the std of the (raw image - fisrt order smoothing)
%   OUTPUTS: 
%       locations:      list of the local maximums [y_top, x_left, frame#]
%       subregion:      sub images cropped around each local maximum
%       subvar:         sub var map cropped around each local maximum
%       ims:            uint16, unprocessed images due to the RAM limit
%       endfrm:         index of the last frame w.r.t the raw image stack after calling this function.
% 




%% ============================ Initialization ============================ 
% Cap of number of PSFs reserved on RAM
Cap = floor(buffer / 4 / (boxsz * boxsz) / 2); 

% 'rafius' of the box
if mod(boxsz, 2) == 0
    error('boxsz must be odd integer.');
end
boxr = (boxsz - 1) / 2; 
lclMax_mask = true(boxsz, boxsz);
lclMax_mask(boxr + 1, boxr + 1) = false;
[imszy, imszx, numImg] = size(ims);


%% ================== 1/var weighted normalization factor =================
switch filtmethod
    case 'box'
        norm1 = imboxfilt(1./var, [2*round(PSFsigma)+1, 2*round(PSFsigma)+1]);
        norm2 = imboxfilt(1./var, [4*round(PSFsigma)+1, 4*round(PSFsigma)+1]);
    case 'gauss'
        norm1 = imgaussfilt(1./var, PSFsigma);
        norm2 = imgaussfilt(1./var, 2*PSFsigma);
    case 'wvlet'
        K = getBSplineKernel(boxsz, 2, 3, 2*PSFsigma);
        norm1 = conv2(1./var, K(1).kernel, 'same');
        norm1 = conv2(norm1, K(1).kernel', 'same');
        norm2 = conv2(1./var, K(2).kernel, 'same');
        norm2 = conv2(norm2, K(2).kernel', 'same');
    otherwise
        error('please choose only from "box", "gauss", or "wvlet"');
end

%% ============= smooth, threshold and segment frame-by-frame =============
locations = zeros(Cap, 3, 'uint16');
subregion = zeros(boxsz, boxsz, Cap, 'single');
subvar = zeros(boxsz, boxsz, Cap, 'single');

ims_res = [];
numOfPeaks = 0;
endfrm = startfrm;
for ff = 1 : numImg
    
    if numOfPeaks > Cap
        ims_res = ims(:, :, ff : end);
        break;
    end
    
    % Always read the 1st image and empty it after segmentation
    I0 = (double(ims(:, :, ff)) - offset) ./ gain;
    endfrm = endfrm + 1;
    
    % ---------------------- 1/var weighted smoothing ---------------------
    switch filtmethod
        case 'box'
            smoothing1 = imboxfilt(I0./var, [2*round(PSFsigma)+1, 2*round(PSFsigma)+1]);
            smoothing2 = imboxfilt(I0./var, [4*round(PSFsigma)+1, 4*round(PSFsigma)+1]);
        case 'gauss'
            smoothing1 = imgaussfilt(I0./var, PSFsigma);
            smoothing2 = imgaussfilt(I0./var, 2*PSFsigma);
        case 'wvlet'
            smoothing1 = conv2(I0./var, K(1).kernel, 'same');
            smoothing1 = conv2(smoothing1, K(1).kernel', 'same');
            smoothing2 = conv2(I0./var, K(2).kernel, 'same');
            smoothing2 = conv2(smoothing2, K(2).kernel', 'same');
    end
    I1 = smoothing1 ./ norm1;
    I2 = smoothing2 ./ norm2;
    
    % --------------------------- thresholding ----------------------------
    if strcmp(thresh, 'std')
        F1 = I0 - I1;
        thresh = std(F1(:));
    end
    F2 = I1 - I2;
    im_max = (F2 > imdilate(F2, lclMax_mask)) & (F2 > thresh);
    
    % ----------------------- seeking local maximum -----------------------
    idx = find(im_max);
    idx_y = mod((idx - 1), imszy) + 1; 
    idx_x = floor((idx - 1) / imszy) + 1;
    
    % exclude spots near the edges
    atedge = idx_x < boxr + 1 | idx_x > imszx - boxr - 1 | idx_y < boxr + 1 | idx_y > imszy - boxr - 1;
    idx_y(atedge, :) = []; 
    idx_x(atedge, :) = [];
    
    % --------------------- sort into location list -----------------------
    dum = size(idx_y, 1);
    locations(numOfPeaks + 1 : numOfPeaks + dum, 1) = idx_y - boxr - 1; % y_top
    locations(numOfPeaks + 1 : numOfPeaks + dum, 2) = idx_x - boxr - 1; % x_left
    locations(numOfPeaks + 1 : numOfPeaks + dum, 3) = endfrm; % frame number (w.r.t the raw image stack)
    for i = 1 : dum
        tmp_data = I0(idx_y(i)-boxr : idx_y(i)+boxr, idx_x(i)-boxr : idx_x(i)+boxr);
        tmp_var = var(idx_y(i)-boxr : idx_y(i)+boxr, idx_x(i)-boxr : idx_x(i)+boxr);
        tmp_gain = gain(idx_y(i)-boxr : idx_y(i)+boxr, idx_x(i)-boxr : idx_x(i)+boxr);
        subregion(:, :, numOfPeaks + i) = single(tmp_data);
        subvar(:, :, numOfPeaks + i) = single(tmp_var ./ tmp_gain ./ tmp_gain);
    end
    numOfPeaks = numOfPeaks + dum;
 
end

%% ======================= Release Extra Space ============================
if numOfPeaks < Cap
    locations(numOfPeaks + 1 : Cap, :) = [];
    subregion(:, :, numOfPeaks + 1 : Cap) = [];
    subvar(:, :, numOfPeaks + 1 : Cap) = [];
end
