clearvars
clc
fclose('all');
addpath('C:\Users\ACLSETPC009\Documents\MATLAB\z_STORM_Analysis_052021CL\MATLAB\Correlation\smPairCorr\CoordinatesBased');
addpath('C:\Users\ACLSETPC009\Documents\MATLAB\z_STORM_Analysis_052021CL\MATLAB\Correlation\ROIRender');

ColorMode = 'RG'; % the autocorrelation will be calculated for all the colors specified here;
                    % the crosscorrelation will be calculated betweeen every two of the specified color here
                    % Only auto-correlation will be calculated if there is only one color specified
                    
Camera_pxsz = 65.39; % nm 
Target_pxsz = 5; % nm

precision_bin = 0:0.05:0.4*Camera_pxsz; % presicion bigger than 0.2*original_pxsz are already filtered out
CorrLength = 64;%normally this is at 64
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

uncertsumm = struct([]);
for ii = 1 : size(ColorMode, 2)
    uncertsumm(ii).fid = fopen([path0 '\UncertSumm_' cl(ii).lname '.uncert'],'w');
end

isSampleDir = [SampleList(:).isdir];
SampleNames = {SampleList(isSampleDir).name}';
SampleNames(ismember(SampleNames, {'.','..','map','Analysis','preview'})) = [];

for s = 1 : numel(SampleNames)
    
    Precision = zeros(size(precision_bin, 2), size(ColorMode, 2));
    PrecisionFit = zeros(size(precision_bin, 2), size(ColorMode, 2));
    
    for ii = 1 : size(ColorMode, 2)
        if ~exist([path0 '\' SampleNames{s} '\PairCorrelation\' cl(ii).lname],'dir')
            mkdir([path0 '\' SampleNames{s} '\PairCorrelation\' cl(ii).lname]);
        end
    end
    
    for ii = 1 : size(ColorMode, 2)
        for jj = ii+1 : size(ColorMode, 2)
            if ~exist([path0 '\' SampleNames{s} '\PairCorrelation\' cl(ii).sname 'and' cl(jj).sname],'dir')
                mkdir([path0 '\' SampleNames{s} '\PairCorrelation\' cl(ii).sname 'and' cl(jj).sname]);
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
            [roi, PrecisionCnts] = ROIRender_v22112020(reports, vnRectBounds(r,:), precision_bin, Camera_pxsz, Target_pxsz);
            
            % ==================== ROI AutoCorrelation ====================
            for ii = 1 : size(ColorMode, 2)
                
                Precision(:, ii) = Precision(:, ii) + PrecisionCnts(ii).cnts';
                
                CorrMat = struct([]);
                
                [bins, ~, Corrtmp] = smPairCorrelation_CPP([roi(ii).x, roi(ii).y], [roi(ii).x, roi(ii).y], roi(ii).roiszx, roi(ii).roiszy);
                clear mex
                ind_NaN = isnan(Corrtmp);
                Corrtmp(ind_NaN, :) = [];
                
                if ~isempty(Corrtmp)
                    bins(ind_NaN, :) = [];
                    Corr = [bins, Corrtmp];
                else
                    Corr = [bins, ones(size(bins, 1), 1)];
                end
                
                CorrMat(1).density_1 = roi(ii).density;
                CorrMat(1).density_2 = roi(ii).density;
                CorrMat(1).Corr = Corr;
                
                save([path0 '\' SampleNames{s} '\PairCorrelation\' cl(ii).lname '\' SpoolName '_' cl(ii).lname '_' num2str(r) '.mat'],'CorrMat');
                clear CorrMat
                
            end
            
            % ===================== ROI CrossCorrelation ==================
            for ii = 1 : size(ColorMode, 2)-1
                for jj = ii+1 : size(ColorMode, 2)
                    CorrMat = struct([]);
                    
                    [bins, ~, Corrtmp] = smPairCorrelation_CPP([roi(ii).x, roi(ii).y], [roi(jj).x, roi(jj).y], roi(ii).roiszx, roi(ii).roiszy);
                    clear mex
                    ind_NaN = isnan(Corrtmp);
                    Corrtmp(ind_NaN, :) = [];
                
                    if ~isempty(Corrtmp)
                        bins(ind_NaN, :) = [];
                        Corr = [bins, Corrtmp];
                    else
                        Corr = [bins, ones(size(bins, 1), 1)];
                    end
                    
                    CorrMat(1).density_1 = roi(ii).density;
                    CorrMat(1).density_2 = roi(jj).density;
                    CorrMat(1).Corr = Corr;
                    save([path0 '\' SampleNames{s} '\PairCorrelation\' cl(ii).sname 'and' cl(jj).sname '\' SpoolName '_' cl(ii).sname cl(jj).sname '_' num2str(r) '.mat'],'CorrMat');
                    clear CorrMat
                end
            end
            
            toc;
        end
    end
    
    options = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','MaxFunEvals',40000,'MaxIter',2000,'TolFun',1e-8,'TolX',1e-8);
    
    for ii = 1 : size(ColorMode, 2)
        indtmp = find(Precision(:, ii)>0, 1, 'last');
        XDATA = precision_bin(1 : indtmp);
        XDATA = XDATA';
        YDATA = Precision(1:indtmp, ii);
        
        [~, tmp_I] = max(YDATA);
        % skewgauss fit of precision peak, % x0 = [a, xc, sigma, alpha, y0]
        para0 = [(max(YDATA) - min(YDATA))*sqrt(2*pi)*10, tmp_I*0.05, 10, 3, min(YDATA)];
        parafit = lsqcurvefit(@skewguassfit, para0, XDATA, YDATA, [], [], options);
        if parafit(2) <= 0 || parafit(3) <= 0
            parafit(2) = sum(XDATA.*YDATA)/sum(YDATA);
            parafit(3) = sqrt(sum((XDATA.*YDATA - parafit(2)).^2)/(sum(YDATA)-1));
        end
        fittedDev = skewguassfit(parafit, XDATA); 
        PrecisionFit(1:indtmp, ii) = fittedDev; 
        
        fprintf(uncertsumm(ii).fid, '%s %f %f %f\r\n', SampleNames{s}, parafit(2), parafit(3), parafit(4));
        
        h = 1;
        figure(h)
        plot(XDATA, YDATA, 'ko'); title([cl(ii).lname ' Dev-Dev hist']);
        hold on
        plot(XDATA, fittedDev,'r-','LineWidth',2);
        hold off
        
        saveas(h, [path0 '\' SampleNames{s} '\PairCorrelation\' cl(ii).lname '_Uncert_hist.jpg'],'jpg');
        clf(h);    
    end
    
    Uncerthist = [precision_bin' Precision, PrecisionFit];
    save([path0 '\' SampleNames{s} '\PairCorrelation\Uncert.hist'],'Uncerthist','-ASCII');
end

for ii = 1 : size(ColorMode, 2)
    fclose(uncertsumm(ii).fid);
end
