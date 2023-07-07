function [roi, PrecisionCnts] = ROIRender_v1129(reports, vnRectBounds, precision_bin, original_pxsz, target_pxsz, mode)

% ================================ INPUT ==================================
% reports = sCMOS MLE results table, each element contains result table of one color channel with col1 = x, col2 = y. Unit: original_pxsz
% vnRectBounds = [top left bottom right]. Unit: original_pxsz
% precision_bin = pre-assigned bins for precision statistics. Unit: nm
% original_pxsz = camera_pxsz. Unit: nm
% target_pxsz = the pxsz for the output roi. Unit: nm
% mode = 'MFA' or 'SFA'. It affects where to read the precision information from the reports table
%
% ================================ OUTPUT =================================
% roi: structure with elements representing different color channels, each of which contains coordinates (nm), roisz (nm), density (cnts/nm^2), and rendered image (target_pxsz)
% PrecisionCnts = each element represents histogrammed precisions of one color channel
%
% =========================== Required Functions ==========================
% grouping_result = overctgroup([x_ini, y_ini], precision.^2, frm_ini); % Grouping artificial blinkings
%
%
% ======================== Calculating Roiszs =============================
vnRectBounds = double(vnRectBounds);
precision_bin = [precision_bin precision_bin(end)+precision_bin(end)-precision_bin(end-1)];
roiszx = vnRectBounds(4) - vnRectBounds(2); % Unit: original_pxsz
roiszy = vnRectBounds(3) - vnRectBounds(1); % Unit: original_pxsz

% =================== Reading Coordinates and Precisions ==================
switch mode
    case 'MFA'
        precision_col = 4;
    case 'SFA'
        precision_col = 5;
    otherwise
        error('reconstruction method (SFA or MFA) is not identified');
end

Source = struct([]);
for ii = 1 : size(reports, 2)
    Source(ii).x = reports(ii).table(:,1); % Unit: original_pxsz
    Source(ii).y = reports(ii).table(:,2); % Unit: original_pxsz
    Source(ii).precision = reports(ii).table(:, precision_col); % Unit: original_pxsz
    Source(ii).frm = reports(ii).table(:,8);
end

% ======== Masking coordinates and precisions that within the roi =========
roi = struct([]);
PrecisionCnts = struct([]);

for ii = 1 : size(reports, 2) % ii indexing color channels
     
    Mask = Source(ii).x>vnRectBounds(2) & Source(ii).x<vnRectBounds(4) & Source(ii).y>vnRectBounds(1) & Source(ii).y<vnRectBounds(3);
    
    x_ini = Source(ii).x(Mask) - vnRectBounds(2); % Unit: original_pxsz
    y_ini = Source(ii).y(Mask) - vnRectBounds(1); % Unit: original_pxsz
    precision_ini = Source(ii).precision(Mask); % Unit original_pxsz
    frm_ini = Source(ii).frm(Mask);
    
    precision_roi = precision_ini.*original_pxsz; % Unit: nm
    PrecisionCnts(ii).cnts = histcounts(precision_roi, precision_bin); % precision statistics

    grouping_result = overctgroup([x_ini, y_ini], precision_ini.^2, frm_ini); % Grouping artificial blinkings
    grouping_result(isnan(grouping_result(:, 1)) | isnan(grouping_result(:, 2)), :) = [];
    roi(ii).x = grouping_result(:, 1) .* original_pxsz; % Unit: nm
    roi(ii).y = grouping_result(:, 2) .* original_pxsz; % Unit: nm
    
    img = zeros(ceil(roiszy*original_pxsz/target_pxsz), ceil(roiszx*original_pxsz/target_pxsz), 'uint8');
    indx = floor(grouping_result(:, 1).*original_pxsz./target_pxsz);
    indy = floor(grouping_result(:, 2).*original_pxsz./target_pxsz);
    idx = indx.*ceil(roiszy*original_pxsz/target_pxsz) + indy + 1;
    img(idx) = 1;
    roi(ii).img = img; % Unit: target_pxsz
    
    roi(ii).roiszx = roiszx * original_pxsz; % Unit: nm
    roi(ii).roiszy = roiszy * original_pxsz; % Unit: nm
    roi(ii).density = size(grouping_result(:, 1), 1) / (roiszx*roiszy*original_pxsz^2); % Unit: cnts/nm^2
    
end