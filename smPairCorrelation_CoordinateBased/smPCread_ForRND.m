clearvars
clc
fclose('all');
addpath('C:\Users\ACLSETPC009\Documents\MATLAB\z_STORM_Analysis_052021CL\MATLAB\Correlation\smPairCorr\CoordinatesBased');
addpath('C:\Users\ACLSETPC009\Documents\MATLAB\z_STORM_Analysis_052021CL\MATLAB\Correlation\ROIRender');
begin = 2;
ColorMode = 'DBRG';
keyword = 'spool_';

options_CC = optimoptions(@lsqcurvefit, 'Algorithm','trust-region-reflective','MaxFunEvals',40000, 'MaxIter',2000,'TolFun',1e-8,'TolX',1e-8);
% -------------------------------------------------------------------------

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
    
    CrossLb = [0, 0, 0, 0.95];
    CrossUb = [inf, inf, inf, inf];
    for ii = 1 : size(ColorMode, 2)
        for jj = ii+1 : size(ColorMode, 2)
            
            pth = [path0 '\' SampleNames{s} '\PairCorrelation_ForRND\' cl(ii).sname 'and' cl(jj).sname '\'];
            CrossFit = zeros(1, 3);
        
            fid_fit = fopen([path0 '\' SampleNames{s} '\PairCorrelation_ForRND\Cross' cl(ii).sname cl(jj).sname 'fit.txt'],'w');
            fprintf(fid_fit,'%-20s %-12s %-12s %-12s %-12s\r\n', 'ROI', 'Autorho_1', 'Autorho_2', 'Amp_norm', 'Amp_unnorm');
            fitFormat = '%-20s %-12.4E %-12.4E %-12.4E %-12.4E\r\n';
            
            fid_fitpara = fopen([path0 '\' SampleNames{s} '\PairCorrelation_ForRND\Cross' cl(ii).sname cl(jj).sname 'fitpara.txt'],'w');
            fprintf(fid_fitpara,'%-20s %-12s %-12s %-12s %-12s\r\n', 'ROI', 'Amp', 'xc', 'sigma', 'y0');
            fitparaFormat = '%-20s %-12.4E %-12.4E %-12.4E %-12.4E\r\n';
            
            fid_fitparaSD = fopen([path0 '\' SampleNames{s} '\PairCorrelation_ForRND\Cross' cl(ii).sname cl(jj).sname 'fitparaSD.txt'],'w');
            fprintf(fid_fitparaSD,'%-20s %-12s %-12s %-12s %-12s %-12s\r\n', 'ROI', 'Amp', 'xc', 'sigma', 'y0');
            fitparaSDFormat = '%-20s %-12.4E %-12.4E %-12.4E %-12.4E\r\n';
        
            CorrFile = dir([pth '*.mat']);
        
            for sp = 1 : numel(CorrFile)
        
                fname = CorrFile(sp).name;
                indtmp1 = strfind(fname, keyword);
                spoolnum = sscanf(fname(indtmp1+length(keyword) : end), '%d', 1);
                indtmp2 = strfind(fname, [cl(ii).sname cl(jj).sname]);
                roinum = sscanf(fname(indtmp2+length(cl(ii).sname)+length(cl(jj).sname)+1 : end), '%d', 1);
            
                roiname = [keyword '_' num2str(spoolnum) '_' num2str(roinum)];
            
                CorrMatFile = load([pth '\' CorrFile(sp).name]);
                Correlation = double(CorrMatFile.CorrMat.Corr);
                density_1 = double(CorrMatFile.CorrMat.density_1);
                density_2 = double(CorrMatFile.CorrMat.density_2);
                Correlation(isnan(Correlation(:,2)), :) = [];
                Correlation(isinf(Correlation(:,2)), :) = [];
            
                if size(Correlation, 1) > 0
                    para0 = [max(Correlation(:,2)), 100, 100, 1];
                    XDATA = Correlation(:, 1);
                    YDATA = Correlation(:, 2);
                    [parafit,~,Res,~,~,~,J] = lsqcurvefit(@crossfit, para0, XDATA, YDATA, CrossLb, CrossUb, options_CC);
                    YDATA_fit = crossfit(parafit, XDATA);
                    ChiSquare = sum((YDATA - YDATA_fit).^2./YDATA_fit);
                    pval = 1 - chi2cdf(ChiSquare, size(YDATA, 1)-4);
                
                    ci = nlparci(parafit, Res, 'jacobian', J);
                    fprintf(fid_fitpara, fitparaFormat, roiname, parafit);
                    fprintf(fid_fitparaSD, fitparaSDFormat, roiname, (ci(:,2) - ci(:,1))./2');
                
                    v = exp(-0.5*parafit(2)^2/parafit(3)^2);
                    corrl = parafit(3)*v*sqrt(2/pi) + parafit(2)*sqrt(1-v);
                    probBG = crossfit(parafit, corrl);
                    
                    CrossFit(1) = parafit(1);
                    CrossFit(2) = corrl;
                    CrossFit(3) = pval;
                    
                    fprintf(fid_fit, fitFormat, roiname, density_1, density_2, CrossFit(1), density_1*density_2*CrossFit(1));
                    
                    h = 1;
                    figure(h);
                    semilogx(XDATA, YDATA, 'ko');
                    hold on
                    semilogx(XDATA, YDATA_fit,'r-','LineWidth',1.5);
                    hold off
                    saveas(h, [pth '\' fname(1:end-5) num2str(roinum) '.jpg']);
                    clf(h);
                    
                    results = [Correlation(:, 1:2), YDATA_fit];
                    save([pth '\' fname(1:end-5) num2str(roinum) '_fit.txt'],'results','-ASCII');
                end
            end
            fclose(fid_fit);
            fclose(fid_fitpara);
            fclose(fid_fitparaSD);
        end
    end
end

function [calcu_ydata] = crossfit(para, xdata)
% expression : y = y0 + a*exp(-((x-c)/r)^2/2)+ a*exp(-((x+c)/r)^2/2)
% para = [a c r y0] for cross fit
calcu_ydata = para(4) + para(1).*exp(-0.5.*(xdata-para(2)).^2./para(3)^2) + para(1).*exp(-0.5.*(xdata+para(2)).^2./para(3)^2);
end