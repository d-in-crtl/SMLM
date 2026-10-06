function [roi_coors] = getROIcoors(reports, roiObj, Camera_pxsz, Image_pxsz, ReconMode)

%{ 
 ================================= INPUT ==================================
 
 reports = sCMOS MLE results table, each element contains result table of one color channel with col1 = x, col2 = y. Unit: Camera_pxsz
 roiObj = an object that contains information of roi. Unit: Image_pxsz
 Camera_pxsz, Unit: nm
 Image_pxsz, Unit: nm (This is the image where the roi is obtained)
 mode = 'MFA' or 'SFA'. It affects where to read the precision information from the reports table

 ================================ OUTPUT =================================
 
 coordinates within the roi

 =========================== Required Functions ==========================
 
 grouping_result = overctgroup([x_ini, y_ini], precision.^2, frm_ini); % Grouping artificial blinkings
 roiObj = readroi(zipFileName) % reading roi information

%}

% =================== Reading Coordinates and Precisions ==================
switch ReconMode
    case 'MFA'
        precision_col = 4;
    case 'SFA'
        precision_col = 5;
    case 'MFA_new'
        precision_col = 9;        
    otherwise
        error('reconstruction method (SFA or MFA) is not identified');
end

% ======== Masking coordinates and precisions that within the roi =========
roi_coors = struct([]);

for ii = 1 : numel(roiObj) % ii indexing color channels
    
    % Note to match the unit of coordinates in reports (Camera_pxsz) and that in the roi measures (Image_pxsz)   
    inside = inpolygon(reports(:, 1).*Camera_pxsz, reports(:, 2).*Camera_pxsz, roiObj(ii).bdry(:, 1).*Image_pxsz, roiObj(ii).bdry(:, 2).*Image_pxsz);
    
    x_ini = reports(inside, 1) .* Camera_pxsz; % Unit nm
    y_ini = reports(inside, 2) .* Camera_pxsz; % Unit nm
    precision_ini = reports(inside,  precision_col) .* Camera_pxsz; % Unit nm
    frm_ini = reports(inside, 8);

    grouping_result = overctgroup([x_ini, y_ini], precision_ini.^2, frm_ini); % Grouping artificial blinkings
    roi_coors(ii).x = grouping_result(:, 1); % Unit: nm
    roi_coors(ii).y = grouping_result(:, 2); % Unit: nm
   
    roi_coors(ii).area = roiObj(ii).area * Image_pxsz * Image_pxsz ; % Unit: nm^2
    roi_coors(ii).density = size(grouping_result(:, 1), 1) / roi_coors(ii).area; % Unit: cnts/nm^2
    
end