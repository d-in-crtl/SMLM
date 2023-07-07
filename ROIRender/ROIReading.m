function [vnRectBounds] = ROIReading(zipFileName)

    [pth, tempdir, ~] = fileparts(zipFileName);
    tmpROIdir = [pth '\' tempdir];

    if exist(tmpROIdir,'dir') == 7
        rmdir(tmpROIdir, 's');
        mkdir(tmpROIdir);
    else
        mkdir(tmpROIdir);
    end
    unzip(zipFileName, tmpROIdir);

    ROIs = dir([tmpROIdir '\*.roi']);
    vnRectBounds = zeros(numel(ROIs), 4, 'int16');

    for i = 1 : numel(ROIs)
    
        ROIname = [tmpROIdir '\' ROIs(i).name];
        fidROI = fopen(ROIname, 'r', 'ieee-be');
    
        strMagic = fread(fidROI, [1 4], '*char'); % read the strMagic
        if (~isequal(strMagic, 'Iout'))
            error('ReadImageJROI:FormatError', '*** ReadImageJROI: The file was not an ImageJ ROI format.');
        end
        Version = fread(fidROI, 1, 'int16'); % read version
        TypeID = fread(fidROI, 1, 'uint8'); % read ROI type
        fseek(fidROI, 1, 'cof'); % skip a byte
        vnRectBounds(i, :) = fread(fidROI, [1 4], 'int16'); % read rectangular bounds [top, left, bottom, right]
        fclose(fidROI);
        
    end
    
    if max(vnRectBounds(:)) > 1200 || min(vnRectBounds(:)) < 0
        error('ROI trying to exceed the image size');
    end
    
    rmdir(tmpROIdir, 's');
end


