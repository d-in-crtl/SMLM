clearvars
clc
fclose('all');
addpath('C:\Users\ACLSETPC009\Documents\MATLAB\z_STORM_Analysis_052021CL\MATLAB\Correlation\smPairCorr\CoordinatesBased');
addpath('C:\Users\ACLSETPC009\Documents\MATLAB\z_STORM_Analysis_052021CL\MATLAB\Correlation\ROIRender');

ColorMode = 'DBRG'; % the autocorrelation will be calculated for all the colors specified here;
                    % the crosscorrelation will be calculated betweeen every two of the specified color here
                    % Only auto-correlation will be calculated if there is only one color specified
Camera_pxsz = 65.39; %73.3; % nm
Target_pxsz = 5; % nm

precision_bin = 0:0.05:0.46*Camera_pxsz; % presicion bigger than 0.2*original_pxsz are already filtered out
CorrLength = 64;

kw1 = 'spool_';
kw2 = '_roi_'; 
% =========================================================================

cl = struct([]);
for ii = 1 : size(ColorMode, 2)
    switch ColorMode(ii)
        case 'D'
            cl(ii).lname = 'DarkRed';
            cl(ii).sname = 'D';
        case 'R'
            cl(ii).lname = 'RED';
            cl(ii).sname = 'R';
        case 'G'
            cl(ii).lname = 'GREEN';
            cl(ii).sname = 'G';
        case 'B'
            cl(ii).lname = 'BLUE';
            cl(ii).sname = 'B';
        otherwise
            error('ColorMode is not properly specified: D = A750, R = A647, G = A568, and B = A488');
    end
end

path0 = uigetdir('', 'Choose the Output Directory');
SampleList = dir(path0);
isSampleDir = [SampleList(:).isdir];
SampleNames = {SampleList(isSampleDir).name}';
SampleNames(ismember(SampleNames, {'.','..','map','Analysis','preview'})) = [];

for s = 1 : numel(SampleNames)
    
    for ii = 1 : size(ColorMode, 2)
        for jj = ii+1 : size(ColorMode, 2)
            if ~exist([path0 '\' SampleNames{s} '\PairCorrelation_ForRND\' cl(ii).sname 'and' cl(jj).sname],'dir')
                mkdir([path0 '\' SampleNames{s} '\PairCorrelation_ForRND\' cl(ii).sname 'and' cl(jj).sname]);
            end
        end
    end
    
    fidlog = fopen([path0 '\' SampleNames{s} '\PairCorrelation_ForRND\randomization.log'], 'w');
    
    AutoFitResults = struct([]);
    for ii = 1 : size(ColorMode, 2)
        fid = fopen([path0 '\' SampleNames{s} '\Auto' cl(ii).lname 'fit.txt']);
        dummy = textscan(fid, '%s %f %f %f %f %f %f %f %f', 'HeaderLines', 1);
        rois = dummy{1};
        AutoFitResults(ii).table = zeros(numel(rois), 3);
        for kk = 1 : numel(rois)
            AutoFitResults(ii).table(kk, :) = [sscanf(rois{kk},'spool__%d_%d', [1 2]), dummy{2}(kk)];
        end
    end
        
    SpoolZips = dir([path0 '\' SampleNames{s} '\*_roi.zip']);
    
    for sp = 1 : numel(SpoolZips)
        
        ZipName = SpoolZips(sp).name;
        indtmp1 = strfind(ZipName, kw1);
        spnum = sscanf(ZipName(indtmp1+length(kw1) : end), '%d', 1);
     
        FullZipName = [path0 '\' SampleNames{s} '\' ZipName];
        [vnRectBounds] = ROIReading(FullZipName); % ROIs must have the same size!!!
        
        SpoolName = ZipName(1:end-8);
        reports = struct([]);
        for ii = 1 : size(ColorMode, 2)
            reports(ii).table = dlmread([path0 '\' SampleNames{s} '\' SpoolName '_' cl(ii).lname '.result']);
        end
        
        
        %---- generates a random sp that different from the current sp ----
        sp_rnd = datasample([1:sp-1, sp+1:numel(SpoolZips)], 1);
        %------------------------------------------------------------------
        ZipName_rnd = SpoolZips(sp_rnd).name;
        indtmp1 = strfind(ZipName_rnd, kw1);
        spnum_rnd = sscanf(ZipName_rnd(indtmp1+length(kw1) : end), '%d', 1);
        
        FullZipName_rnd = [path0 '\' SampleNames{s} '\' ZipName_rnd];
        [vnRectBounds_rnd] = ROIReading(FullZipName_rnd); % ROIs must have the same size!!!
        
        SpoolName_rnd = ZipName_rnd(1:end-8);
        reports_rnd = struct([]);
        for ii = 1 : size(ColorMode, 2)
            reports_rnd(ii).table = dlmread([path0 '\' SampleNames{s} '\' SpoolName_rnd '_' cl(ii).lname '.result']);
        end
        
        
        for r = 1 : size(vnRectBounds,1)
            
            tic; 
            [roi, ~] = ROIRender_v22112020(reports, vnRectBounds(r,:), precision_bin, Camera_pxsz, Target_pxsz);
            
            %---------------------- generates a random r ------------------
            r_rnd = datasample((1:size(vnRectBounds_rnd,1)), 1);
            %--------------------------------------------------------------
            [roi_rnd, ~] = ROIRender_v22112020(reports_rnd, vnRectBounds_rnd(r_rnd,:), precision_bin, Camera_pxsz, Target_pxsz);
            
            msg = [SpoolName '-roi-' num2str(r) ' vs ' SpoolName_rnd '-roi-' num2str(r_rnd)];
            fprintf(fidlog, '%s\r\n', msg);
            disp([num2str(s) ': ' SampleNames{s} '-' msg]);
            
            
            % ===================== ROI CrossCorrelation ==================
            for ii = 1 : size(ColorMode, 2)-1
                for jj = ii+1 : size(ColorMode, 2)
                    ind_i = find(AutoFitResults(ii).table(:,1) == spnum & AutoFitResults(ii).table(:,2) == r);
                    ind_j = find(AutoFitResults(jj).table(:,1) == spnum_rnd & AutoFitResults(jj).table(:,2) == r_rnd);
                    [bins, ~, Corrtmp] = smPairCorrelation_CPP([roi(ii).x, roi(ii).y], [roi_rnd(jj).x, roi_rnd(jj).y], roi(ii).roiszx, roi(ii).roiszy);
                    clear mex;
                    ind_NaN = isnan(Corrtmp);
                    Corrtmp(ind_NaN, :) = [];
                    if ~isempty(Corrtmp)
                        bins(ind_NaN, :) = [];
                        Corr = [bins, Corrtmp];
                    else
                        Corr = [bins, ones(size(bins, 1), 1)];
                    end
                    CorrMat(1).density_1 = AutoFitResults(ii).table(ind_i, 3);
                    CorrMat(1).density_2 = AutoFitResults(jj).table(ind_j, 3);
                    CorrMat(1).Corr = Corr;
                    fprintf(fidlog, '%s %s %f %f\r\n', cl(ii).sname, cl(jj).sname, CorrMat(1).density_1, CorrMat(1).density_2);
                    %disp([CorrMat(1).density_1, CorrMat(1).density_2]);
                    save([path0 '\' SampleNames{s} '\PairCorrelation_ForRND\' cl(ii).sname 'and' cl(jj).sname '\' SpoolName '_' cl(ii).sname cl(jj).sname '_' num2str(r) '.mat'],'CorrMat');
                    clear CorrMat;
                end
            end
            
            clear roi roi_rnd
            toc;
        end
        clear reports reports_rnd
    end
    fclose(fidlog);
end
