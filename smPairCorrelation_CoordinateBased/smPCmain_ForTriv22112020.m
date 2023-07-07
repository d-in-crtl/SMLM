clearvars
clc
fclose('all');
addpath('C:\Users\ACLSETPC009\Documents\MATLAB\z_STORM_Analysis_052021CL\MATLAB\Correlation\smPairCorr\CoordinatesBased');
addpath('C:\Users\ACLSETPC009\Documents\MATLAB\z_STORM_Analysis_052021CL\MATLAB\Correlation\ROIRender');

ColorMode = 'DBRG'; % the autocorrelation will be calculated for all the colors specified here;
                    % the crosscorrelation will be calculated betweeen every two of the specified color here
                    % Only auto-correlation will be calculated if there is only one color specified
                    
Camera_pxsz = 65.39; % nm  
Target_pxsz = 5; % nm

precision_bin = 0:0.05:0.4*Camera_pxsz; % presicion bigger than 0.2*original_pxsz are already filtered out
CorrLength = 200;
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
            if ~exist([path0 '\' SampleNames{s} '\PairCorrelation_ForTri\' cl(ii).sname 'and' cl(jj).sname],'dir')
                mkdir([path0 '\' SampleNames{s} '\PairCorrelation_ForTri\' cl(ii).sname 'and' cl(jj).sname]);
            end
        end
    end
    
    SpoolZips = dir([path0 '\' SampleNames{s} '\*_roi.zip']);
    
    for sp = 1 : numel(SpoolZips)
        
        ZipName = SpoolZips(sp).name;
        
        SpoolName = ZipName(1:end-8);
        reports = struct([]);
        for ii = 1 : size(ColorMode, 2)
            reports(ii).table = dlmread([path0 '\' SampleNames{s} '\' SpoolName '_' cl(ii).lname '.result']);
        end
        
        FullZipName = [path0 '\' SampleNames{s} '\' ZipName];
        [vnRectBounds] = ROIReading(FullZipName);
        
        for r = 1 : size(vnRectBounds,1)
            disp([num2str(s) ': ' SampleNames{s} '-' SpoolName '-roi-' num2str(r)]);
            tic; 
            [roi, ~] = ROIRender_v22112020(reports, vnRectBounds(r,:), precision_bin, Camera_pxsz, Target_pxsz);
            
            % ===================== ROI CrossCorrelation ==================
            for ii = 1 : size(ColorMode, 2)-1
                for jj = ii+1 : size(ColorMode, 2)
                    [bins, ~, Corrtmp] = smPairCorrelation_ForTri([roi(ii).x, roi(ii).y], [roi(jj).x, roi(jj).y], roi(ii).roiszx, roi(ii).roiszy);
                    ind_NaN = isnan(Corrtmp);
                    Corrtmp(ind_NaN, :) = [];
                    if ~isempty(Corrtmp)
                        bins(ind_NaN, :) = [];
                        Corr = [bins, Corrtmp];
                    else
                        Corr = [bins, ones(size(bins, 1), 1)];
                    end
                    save([path0 '\' SampleNames{s} '\PairCorrelation_ForTri\' cl(ii).sname 'and' cl(jj).sname '\' SpoolName '_' cl(ii).sname cl(jj).sname '_' num2str(r) '.ForTri'],'Corr','-ASCII');
                end
            end
            
            clear roi
            toc;
        end
        clear reports
    end
    
end
