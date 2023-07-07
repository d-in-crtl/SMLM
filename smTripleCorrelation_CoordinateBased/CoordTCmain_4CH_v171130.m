clearvars
clc
fclose('all');

addpath('C:\Users\ACLSETPC009\Documents\MATLAB\z_STORM_Analysis_052021CL\MATLAB\Correlation\ROIRender');
addpath('C:\Users\ACLSETPC009\Documents\MATLAB\z_STORM_Analysis_052021CL\MATLAB\Correlation\smTriCorr\CoordinatesBased');

CalMode = 'DBRG';

original_pxsz = 65.39; % nm
rho_res = 5; % nm
MaxRhoSize = 100; % samplesize

mode = 'MFA_new';
clustering = 'N';
grids = 8;

zm = original_pxsz / rho_res;

rho_center = ((1:MaxRhoSize)' - 0.5).*rho_res; 
rho_edge = (0:MaxRhoSize)'.*rho_res;

precision_bin = 0:0.05:0.5*original_pxsz; 

rad_res = pi/45; 
MaxRadSize = ceil(2*pi/rad_res);
rad_res = 2*pi/MaxRadSize;
%rad_center = ((1:MaxRadSize)' - 0.5).*rad_res;


% Read Data

path0 = 'Z:\homes\guptad07\Dipika\062121\Output_2DGauss_MFA_ABFGSB_1';
SampleList = dir(path0);

isSampleDir = [SampleList(:).isdir];
SampleNames = {SampleList(isSampleDir).name}';
SampleNames(ismember(SampleNames, {'.','..','map','Analysis','preview'})) = [];

TripleCorrelation = zeros(MaxRhoSize, MaxRadSize, MaxRhoSize);

for i = 1 : numel(SampleNames)
        
    spools = dir([path0 '/' SampleNames{i} '/*_roi.zip']);
    
    if ~exist([path0 '/' SampleNames{i} '/TripleResults_DBR'],'dir')
        mkdir([path0 '/' SampleNames{i} '/TripleResults_DBR']);
    end
    
    if ~exist([path0 '/' SampleNames{i} '/TripleResults_DBG'],'dir')
        mkdir([path0 '/' SampleNames{i} '/TripleResults_DBG']);
    end
    
    if ~exist([path0 '/' SampleNames{i} '/TripleResults_DRG'],'dir')
        mkdir([path0 '/' SampleNames{i} '/TripleResults_DRG']);
    end
       
    if ~exist([path0 '/' SampleNames{i} '/TripleResults_BRG'],'dir')
        mkdir([path0 '/' SampleNames{i} '/TripleResults_BRG']);
    end
               
               
    for j = 1 : numel(spools)
        
	roi_name = spools(j).name;
        roi_zip = [path0 '/' SampleNames{i} '/' roi_name];
        
        [vnRectBounds] = ROIReading(roi_zip);
        
        fname = roi_name(1:end-8);
        filename_D = [path0 '/' SampleNames{i} '/' fname '_DarkRed.result'];
        filename_R = [path0 '/' SampleNames{i} '/' fname '_RED.result'];
        filename_G = [path0 '/' SampleNames{i} '/' fname '_GREEN.result'];
        filename_B = [path0 '/' SampleNames{i} '/' fname '_BLUE.result'];
               
        reports_D = dlmread(filename_D);
        reports_R = dlmread(filename_R);
        reports_G = dlmread(filename_G);
        reports_B = dlmread(filename_B);
            
        for r = 1 : size(vnRectBounds,1)
	       
	    roiszx = double(vnRectBounds(r,4) - vnRectBounds(r,2)) * original_pxsz;
	    roiszy = double(vnRectBounds(r,3) - vnRectBounds(r,1)) * original_pxsz;

            disp([num2str(i) '-' num2str(j) '-' num2str(r)]);

	    [x_D, y_D, ~, ~] = ROIRender_v0726(reports_D, vnRectBounds(r,:), zm, precision_bin, original_pxsz, mode, clustering, grids);
	    [x_R, y_R, ~, ~] = ROIRender_v0726(reports_R, vnRectBounds(r,:), zm, precision_bin, original_pxsz, mode, clustering, grids);
	    [x_G, y_G, ~, ~] = ROIRender_v0726(reports_G, vnRectBounds(r,:), zm, precision_bin, original_pxsz, mode, clustering, grids);
	    [x_B, y_B, ~, ~] = ROIRender_v0726(reports_B, vnRectBounds(r,:), zm, precision_bin, original_pxsz, mode, clustering, grids);
        
        coor_D = [x_D, y_D] .* original_pxsz;
        coor_R = [x_R, y_R] .* original_pxsz;
        coor_G = [x_G, y_G] .* original_pxsz;
        coor_B = [x_B, y_B] .* original_pxsz;
        
	    if ~isempty(x_R) && ~isempty(x_B) && ~isempty(x_D)
	        tic;
	        [TripleResult_DBR, ~, rad_center] = smTripleCorrelation_v171115_CPU(coor_D, coor_B, coor_R, roiszx, roiszy, MaxRhoSize, rho_res, rad_res);
	        if ~isempty(TripleResult_DBR)
                    save([path0 '/' SampleNames{i} '/TripleResults_DBR/' fname '_triple_' num2str(r) '.mat'],'TripleResult_DBR');
                end
                toc;
            end

	    if ~isempty(x_G) && ~isempty(x_B) && ~isempty(x_D)
	        tic;
	        [TripleResult_DBG, ~, rad_center] = smTripleCorrelation_v171115_CPU(coor_D, coor_B, coor_G, roiszx, roiszy, MaxRhoSize, rho_res, rad_res);
                if ~isempty(TripleResult_DBG)
                    save([path0 '/' SampleNames{i} '/TripleResults_DBG/' fname '_triple_' num2str(r) '.mat'],'TripleResult_DBG');
                end
	        toc;
            end

	    if ~isempty(x_R) && ~isempty(x_G) && ~isempty(x_D)
	        tic;
                [TripleResult_DRG, ~, rad_center] = smTripleCorrelation_v171115_CPU(coor_D, coor_R, coor_G, roiszx, roiszy, MaxRhoSize, rho_res, rad_res);
                if ~isempty(TripleResult_DRG)
                    save([path0 '/' SampleNames{i} '/TripleResults_DRG/' fname '_triple_' num2str(r) '.mat'],'TripleResult_DRG');
                end
                toc;
            end

	    if ~isempty(x_R) && ~isempty(x_G) && ~isempty(x_B)
                tic;
                [TripleResult_BRG, ~, rad_center] = smTripleCorrelation_v171115_CPU(coor_B, coor_R, coor_G, roiszx, roiszy, MaxRhoSize, rho_res, rad_res);
                if ~isempty(TripleResult_BRG)
                    save([path0 '/' SampleNames{i} '/TripleResults_BRG/' fname '_triple_' num2str(r) '.mat'],'TripleResult_BRG');
                end
                toc;
            end
        
        end
    end   
end
