function [subROI, subROIszx, subROIszy] = subROIgroup(roi, roiszx, roiszy, maxCoorsz)
% =============================== INPUT =================================== 
% roi = structure, each element represents one color channel, each element has its own x and y coordinates. Unit: original_pxsz 
% roiszx, roiszy = size of roi. Unit: original_pxsz
% maxCoorsz = number of coordinates in the channel has the most coordinates.

% ============================== OUTPUT ===================================
% subROI: nested structure, with the first lv elements represents divided sub rois, and the second lv elements represent color channels
% each 2nd lv element contains the coordinates. Unit: original_pxsz
% subROIszx and subROIszy: size of subROIs (they are evenly divided). Unit: original_pxsz

subROI = struct([]);

subs = ceil(maxCoorsz / 2E4);

if subs == 1
    subROINum = 1;
    
elseif subs > 1 && subs <= 4
    subROINum = 2;
    
elseif subs > 4 && subs <= 9
    subROINum = 3;
    
else
    subROINum = 0;
end

if subROINum ~= 0
    subROIszx = roiszx / subROINum;
    subROIszy = roiszy / subROINum;
    for col = 0 : subROINum - 1
        for row = 0 : subROINum - 1
            
            subROIIdx = col * subROINum + row + 1;
            
            bnd_l = (col / subROINum) * roiszx;
            bnd_r = ((col + 1) / subROINum) * roiszx;
            bnd_u = (row / subROINum) * roiszy;
            bnd_b = ((row + 1) / subROINum) * roiszy;
            
            for cl = 1 : size(roi, 2)
                Mask = roi(cl).x >= bnd_l & roi(cl).x < bnd_r & roi(cl).y >= bnd_u & roi(cl).y < bnd_b;
                subROI(subROIIdx).CH(cl).x = roi(cl).x(Mask);
                subROI(subROIIdx).CH(cl).y = roi(cl).y(Mask);
                subROI(subROIIdx).CH(cl).counts = size(roi(cl).x(Mask), 1);
            end
            
        end
    end
    
else
    subROIszx = roiszx;
    subROIszy = roiszy;
    for cl = 1 : size(roi, 2)
        subROIIdx = 1;
        subROI(subROIIdx).CH(cl).x = roi(cl).x;
        subROI(subROIIdx).CH(cl).y = roi(cl).y;
        subROI(subROIIdx).CH(cl).counts = size(roi(cl).x, 1);
    end
    
end