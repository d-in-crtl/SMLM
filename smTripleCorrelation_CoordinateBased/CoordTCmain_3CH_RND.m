clearvars
clc
fclose('all');

addpath('C:\Users\rothenberglab\Documents\MATLAB\ROIRender\');
addpath('C:\Users\rothenberglab\Documents\MATLAB\smTriCorr\CoordinatesBased');
addpath('C:\Users\rothenberglab\Documents\MATLAB\smPairCorr\CoordinatesBased');

TriCalMode = 'BRG';
TriDispMode = 'BGR';

CCdim_1 = [TriCalMode(1) TriCalMode(3)];
CCdim_2 = [TriCalMode(1) TriCalMode(2)];
CCdim_3 = [TriCalMode(2) TriCalMode(3)];

switch TriCalMode(1)
    case 'B'
        ACdim_1 = 'BLUE';
    case 'R'
        ACdim_1 = 'RED';
    case 'G'
        ACdim_1 = 'GREEN';
    case 'D'
        ACdim_1 = 'DarkRed';
end
        
switch TriCalMode(2)
    case 'B'
        ACdim_2 = 'BLUE';
    case 'R'
        ACdim_2 = 'RED';
    case 'G'
        ACdim_2 = 'GREEN';
    case 'D'
        ACdim_2 = 'DarkRed';
end

switch TriCalMode(3)
    case 'B'
        ACdim_3 = 'BLUE';
    case 'R'
        ACdim_3 = 'RED';
    case 'G'
        ACdim_3 = 'GREEN';
    case 'D'
        ACdim_3 = 'DarkRed';
end
        
original_pxsz = 65.39; % nm
rho_res = 5; % nm
MaxRhoSize = 100; % samplesize

mode = 'MFA';
clustering = 'N';
grids = 8;

zm = original_pxsz / rho_res;

rho_center = ((1:MaxRhoSize)' - 0.5).*rho_res; 
rho_edge = (0:MaxRhoSize)'.*rho_res;

precision_bin = 0:0.05:0.5*original_pxsz; 

rad_res = pi/45; 
MaxRadSize = ceil(2*pi/rad_res);
rad_res = 2*pi/MaxRadSize;

dim_ini = 80; % nm
dim = round(dim_ini / rho_res);
boxradii = [3, 3, 3];

NormMode = 'const';
PlotMode = 'raw';
normDist = 200;

% Read Data

path0 = 'Z:\homes\yiny02\Data_2019\20200818\Output';
SampleList = dir(path0);

isSampleDir = [SampleList(:).isdir];
SampleNames = {SampleList(isSampleDir).name}';
SampleNames(ismember(SampleNames, {'.','..','map','Analysis','preview'})) = [];

TripleCorrelation = zeros(MaxRhoSize, MaxRadSize, MaxRhoSize);

for i = 1 : numel(SampleNames)
    
    fprintf("%s\n", SampleNames{i});
    if ~exist([path0 '\' SampleNames{i} '\TripleRND_' TriCalMode],'dir')
        mkdir([path0 '\' SampleNames{i} '\TripleRND_' TriCalMode]);
    end
    
    NumofTris = 0;
    fid = fopen([path0 '\' SampleNames{i} '\TripleRND_' TriCalMode '\results_rng' num2str(dim_ini) '.delta'],'w');
    fprintf(fid, '%-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s\r\n', 'ROI', ['Dist_' CCdim_1], ['Dist_' CCdim_2], ['Dist_' CCdim_3], ['delta_' TriCalMode], ['delta_' CCdim_1], ['delta_' CCdim_2], ['delta_' CCdim_3], ['Rho_' TriCalMode(1)], ['Rho_' TriCalMode(2)], ['Rho_' TriCalMode(3)]);
    
    spools = dir([path0 '\' SampleNames{i} '\*_roi.zip']);
               
    for j0 = 1 : numel(spools)
        
        fname = spools(j0).name(1 : end-8);
        
        rndj = datasample([1:j0-1, j0+1:numel(spools)], 2, 'Replace', false);
        j1 = rndj(1);
        j2 = rndj(2);
        
        [vnRectBounds0] = ROIReading([path0 '\' SampleNames{i} '\' spools(j0).name]);
        num_roi0 = size(vnRectBounds0, 1);
        
        [vnRectBounds1] = ROIReading([path0 '\' SampleNames{i} '\' spools(j1).name]);
        num_roi1 = size(vnRectBounds1, 1);
        
        [vnRectBounds2] = ROIReading([path0 '\' SampleNames{i} '\' spools(j2).name]);
        num_roi2 = size(vnRectBounds2, 1);
        
        
        filename_1 = [path0 '\' SampleNames{i} '\' spools(j0).name(1:end-8) '_' ACdim_1 '.result'];
        filename_2 = [path0 '\' SampleNames{i} '\' spools(j1).name(1:end-8) '_' ACdim_2 '.result']; 
        filename_3 = [path0 '\' SampleNames{i} '\' spools(j2).name(1:end-8) '_' ACdim_3 '.result'];
               
        reports_1 = dlmread(filename_1);
        reports_2 = dlmread(filename_2);
        reports_3 = dlmread(filename_3);
            
        for r0 = 1 : num_roi0
            
            r1 = datasample((1 : num_roi1), 1);
            r2 = datasample((1 : num_roi2), 1);
            
            % rois should be of the same size
            roiszx = double(vnRectBounds0(r0,4) - vnRectBounds0(r0,2)) * original_pxsz;
            roiszy = double(vnRectBounds0(r0,3) - vnRectBounds0(r0,1)) * original_pxsz;

            fprintf("Correlating CH1-sp%d-roi%d, CH2-sp%d-roi%d, and CH3-sp%d-roi%d\n", j0, r0, j1, r1, j2, r2);

            [x_1, y_1, ~, ~] = ROIRender_v0726(reports_1, vnRectBounds0(r0,:), zm, precision_bin, original_pxsz, mode, clustering, grids);
            [x_2, y_2, ~, ~] = ROIRender_v0726(reports_2, vnRectBounds1(r1,:), zm, precision_bin, original_pxsz, mode, clustering, grids);
            [x_3, y_3, ~, ~] = ROIRender_v0726(reports_3, vnRectBounds2(r2,:), zm, precision_bin, original_pxsz, mode, clustering, grids);
        
            coor_1 = [x_1, y_1] .* original_pxsz;
            coor_1(isnan(coor_1(:, 1)) | isnan(coor_1(:, 2)), :) = [];
            coor_2 = [x_2, y_2] .* original_pxsz;
            coor_2(isnan(coor_2(:, 1)) | isnan(coor_2(:, 2)), :) = [];
            coor_3 = [x_3, y_3] .* original_pxsz;
            coor_3(isnan(coor_3(:, 1)) | isnan(coor_3(:, 2)), :) = [];

            if ~isempty(x_1) && ~isempty(x_2) && ~isempty(x_3)
                tic;
                [TripleResult, ~, rad_center] = smTripleCorrelation_v171115_CPU(coor_1, coor_2, coor_3, roiszx, roiszy, MaxRhoSize, rho_res, rad_res);
                
                [bin, ~, PairResult_13] = smPairCorrelation_ForTri(coor_1, coor_3, roiszx, roiszy);
                PC_13 = [bin, PairResult_13];
                [~, ~, PairResult_12] = smPairCorrelation_ForTri(coor_1, coor_2, roiszx, roiszy);
                PC_12 = [bin, PairResult_12];
                [~, ~, PairResult_23] = smPairCorrelation_ForTri(coor_2, coor_3, roiszx, roiszy);
                PC_23 = [bin, PairResult_23];
                
                Rho_AC_1 = size(coor_1, 1) / (roiszx * roiszy);
                Rho_AC_2 = size(coor_2, 1) / (roiszx * roiszy);
                Rho_AC_3 = size(coor_3, 1) / (roiszx * roiszy);
                
                [triple_trans, results_tmp, coordinates_tmp, coordinates_tmp_norm] = TripleRead_vDelta(TripleResult, PC_13, PC_12, PC_23, boxradii, dim, rho_res, NormMode, normDist, TriCalMode, TriDispMode);
                for kk = 1 : size(results_tmp, 1)
                    dummy = [results_tmp(kk, :), Rho_AC_1, Rho_AC_2, Rho_AC_3];
                    fprintf(fid, '%-20s %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E\r\n', [num2str(j0) '_' num2str(r0) 'v' num2str(j1) '_' num2str(r1) 'v' num2str(j2) '_' num2str(r2)], dummy);
                end
                
                if ~isempty(TripleResult)
                    save([path0 '\' SampleNames{i} '\TripleRND_' TriCalMode '\' fname '_triple_raw_' num2str(r0) '.mat'],'TripleResult');
                    save([path0 '\' SampleNames{i} '\TripleRND_' TriCalMode '\' fname '_triple_trans_' num2str(r0) '.mat'],'triple_trans');
                end
                toc;
            end
        
        end
    end 
    fclose(fid);
end
