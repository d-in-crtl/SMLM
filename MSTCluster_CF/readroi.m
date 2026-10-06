function [roi] = readroi(zipFileName)
    
    % Note Rotated ROI is not supported
    
    % roi type_ID code
    POLYGON = 0; 
    RECT = 1; 
    ELLIPSE = 2;
    FREEHAND = 7;
    
    % options code
    SUB_PIXEL_RESOLUTION = 128;
    
    % ---------------------------------------------------------------
    [pth, tempdir, ~] = fileparts(zipFileName);
    tmpROIdir = [pth '\' tempdir];

    if exist(tmpROIdir,'dir') == 7
        rmdir(tmpROIdir, 's');
        mkdir(tmpROIdir);
    else
        mkdir(tmpROIdir);
    end
    unzip(zipFileName, tmpROIdir);
    
    % ---------------------------------------------------------------
    
    ROIs = dir([tmpROIdir '\*.roi']);
    roi = struct([]);

    for i = 1 : numel(ROIs)
    
        ROIname = [tmpROIdir '\' ROIs(i).name];
        fidROI = fopen(ROIname, 'r', 'ieee-be');
    
        strMagic = fread(fidROI, [1 4], '*char'); % read the strMagic
        if (~isequal(strMagic, 'Iout'))
            error('ReadImageJROI:FormatError', '*** ReadImageJROI: The file was not an ImageJ ROI format.');
        end
        
        roi(i).Version = fread(fidROI, 1, 'int16'); % read version
        roi(i).TypeID = fread(fidROI, 1, 'uint8'); % read ROI type
        
        if roi(i).TypeID == POLYGON || roi(i).TypeID == FREEHAND || roi(i).TypeID == RECT || roi(i).TypeID == ELLIPSE
            
            fseek(fidROI, 1, 'cof'); % skip a byte
            roi(i).tlbr = fread(fidROI, [1 4], 'int16'); % read rectangular bounds [top, left, bottom, right]
            
            if roi(i).TypeID == POLYGON || roi(i).TypeID == FREEHAND
                
                roi(i).NumPx = fread(fidROI, 1, 'uint16'); % read number of pixels within the roi
                fseek(fidROI, 30, 'cof'); % jump to read subtype
                
                roi(i).subtype = fread(fidROI, 1, 'int16'); % read the subtype
                if roi(i).subtype ~= 0
                    error(['ROI reading error: subtype of roi-' num2str(i) ' is not supported']);
                end
                
                roi(i).options = fread(fidROI, 1, 'int16'); % read the options
                if bitand(roi(i).options, SUB_PIXEL_RESOLUTION)
                    pixel_type = 'single';
                else
                    pixel_type = 'int16';
                end
                fseek(fidROI, 12, 'cof'); % jump to read pixels
                
                PxX = fread(fidROI, [roi(i).NumPx 1], pixel_type);
                PxY = fread(fidROI, [roi(i).NumPx 1], pixel_type);
                [inds, roi(i).area] = boundary(PxX, PxY);
                roi(i).bdry = [PxX(inds) + roi(i).tlbr(2), PxY(inds) + roi(i).tlbr(1)];
                
            elseif roi(i).TypeID == RECT
                roi(i).bdry = [roi(i).tlbr(2), roi(i).tlbr(1);
                               roi(i).tlbr(2), roi(i).tlbr(3);
                               roi(i).tlbr(4), roi(i).tlbr(3);
                               roi(i).tlbr(4), roi(i).tlbr(1);
                               roi(i).tlbr(2), roi(i).tlbr(1)];
                roi(i).area = (roi(i).tlbr(4) - roi(i).tlbr(2)) * (roi(i).tlbr(3) - roi(i).tlbr(1));
                
            elseif roi(i).TypeID == ELLIPSE
                centerX = (roi(i).tlbr(4) + roi(i).tlbr(2)) / 2;
                centerY = (roi(i).tlbr(3) + roi(i).tlbr(1)) / 2;
                armX = (roi(i).tlbr(4) - roi(i).tlbr(2)) / 2;
                armY = (roi(i).tlbr(3) - roi(i).tlbr(1)) / 2;
                
                theta = linspace(0, 2*pi, 37);
                roi(i).bdry = [centerX + armX.*cos(theta'), centerY + armY.*sin(theta')];
                roi(i).area = pi * armX * armY;
                
            end
            
        else
            error('roi type not supported. Current version supports polygon, rectangle, ellipse, and freehand');
        end
        
        fclose(fidROI);
        
    end
    
    rmdir(tmpROIdir, 's');

end


