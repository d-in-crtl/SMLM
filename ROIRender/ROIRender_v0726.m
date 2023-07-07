function [x_out, y_out, roi, precision_roi_cnts] = ROIRender_v0726(reports, vnRectBounds, zm, precision_bin, original_pxsz, mode, clustering, grids)
% vnRectBounds = [top left bottom right]
vnRectBounds = double(vnRectBounds);
precision_bin = [precision_bin precision_bin(end)+precision_bin(end)-precision_bin(end-1)];
roiszx = vnRectBounds(4) - vnRectBounds(2);
roiszy = vnRectBounds(3) - vnRectBounds(1);

roi = zeros(ceil(roiszy*zm), ceil(roiszx*zm), 'uint16');

x = reports(:,1);
y = reports(:,2);
frm = reports(:,8);
switch mode
    case 'MFA'
        precision = reports(:,4);
    case 'SFA'
        precision = reports(:,5);
    case 'MFA_new'
        precision = (reports(:, end-2) + reports(:, end-1)) ./ 2;
    otherwise
        error('reconstruction method (SFA or MFA) is not identified')
end       
mask = x>vnRectBounds(2) & x<vnRectBounds(4) & y>vnRectBounds(1) & y<vnRectBounds(3);
precision_roi = precision(mask).*original_pxsz;
precision_roi_cnts = histcounts(precision_roi, precision_bin);

x_ini = x(mask);
y_ini = y(mask);           
frm_ini = frm(mask);

x_ini = x_ini - vnRectBounds(2);
y_ini = y_ini - vnRectBounds(1); 

grouping_result = overctgroup([x_ini, y_ini], precision.^2, frm_ini);
x_roi = grouping_result(:, 1);
y_roi = grouping_result(:, 2);
var_roi = grouping_result(:, 3);

if clustering == 'Y'
    
    gridxsz = roiszx / grids;
    gridysz = roiszy / grids;
    clustering_result = [0, 0, 0, 0];
    for i = 1 : grids
        for j = 1 : grids
            mask_tmp = x_roi > (j-1)*gridxsz & x_roi <= j*gridxsz & y_roi > (i-1)*gridysz & y_roi <= i*gridysz;
            clustering_result_tmp = FixedRadiusGrouping([x_roi(mask_tmp), y_roi(mask_tmp)], var_roi(mask_tmp), mean(sqrt(var_roi(mask_tmp)), 1));
            clustering_result = cat(1, clustering_result, clustering_result_tmp);
        end
    end
    clustering_result(1, :) = [];
    
    x_out = clustering_result(:, 1);
    y_out = clustering_result(:, 2);
    
else
    x_out = x_roi;
    y_out = y_roi;
end
    
indx = floor(x_out.*zm);
indy = floor(y_out.*zm);

idx = indx.*ceil(roiszy*zm) + indy + 1;
roi(idx) = 1;
