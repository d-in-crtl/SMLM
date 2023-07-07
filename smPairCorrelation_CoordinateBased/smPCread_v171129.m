clearvars
clc
fclose('all');
addpath('C:\Users\chelsea_pc\Documents\MATLAB\smPairCorr\CoordinatesBased');
addpath('C:\Users\chelsea_pc\Documents\MATLAB\ROIRender\');

begin = 2;
ColorMode = 'DBRG';
keyword = 'spool_';

options_AC = optimoptions(@lsqcurvefit, 'Algorithm','trust-region-reflective','MaxFunEvals',40000, 'MaxIter',2000,'TolFun',1e-8,'TolX',1e-8);
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
uncert = struct([]);
for ii = 1 : size(ColorMode, 2)
    uncertfid = fopen([path0 '\UncertSumm_' cl(ii).lname '.uncert'],'r');
    Uncertainty = textscan(uncertfid, '%s %f %f %f');
    uncert(ii).xc = Uncertainty{2};
    uncert(ii).sigma = Uncertainty{3};
end
SampleList = Uncertainty{1};

for s = 1 : numel(SampleList)
    
    for ii = 1 : size(ColorMode, 2)
        
        pth = [path0 '\' SampleList{s} '\PairCorrelation\' cl(ii).lname '\'];
        
        %AutoLb = [0, max(uncert(ii).xc(s)-1.25*uncert(ii).sigma(s),10), 0, min(uncert(ii).xc(s)+1.25*uncert(ii).sigma(s),50), 0.95];
        %AutoUb = [inf, min(uncert(ii).xc(s)+1.25*uncert(ii).sigma(s),50), inf, 300, inf];
        AutoLb = [0, 5, 0, 20, 0.95];
        AutoUb = [inf, 20, inf, 300, inf];
        
        Autofit = zeros(1,7);
        
        fid_fit = fopen([path0 '\' SampleList{s} '\Auto' cl(ii).lname 'fit.txt'],'w');
        fprintf(fid_fit,'%-20s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\r\n', 'ROI', 'raw_rho', 'rho', 'Uncert', 'Radii', 'Compact', 'NumF', 'NumC', 'pval');
        fitFormat = '%-20s %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E\r\n';
            
        fid_fitpara = fopen([path0 '\' SampleList{s} '\Auto' cl(ii).lname 'fitpara.txt'],'w');
        fprintf(fid_fitpara,'%-20s %-12s %-12s %-12s %-12s %-12s\r\n', 'ROI', 'Amp1', 'sigma1', 'Amp2', 'sigma2', 'y0');
        fitparaFormat = '%-20s %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E\r\n';
            
        fid_fitparaSD = fopen([path0 '\' SampleList{s} '\Auto' cl(ii).lname 'fitparaSD.txt'],'w');
        fprintf(fid_fitparaSD,'%-20s %-12s %-12s %-12s %-12s %-12s\r\n', 'ROI', 'Amp1', 'sigma1', 'Amp2', 'sigma2', 'y0');
        fitparaSDFormat = '%-20s %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E\r\n';
        
        CorrFile = dir([pth '*.mat']);
        
        for sp = 1 : numel(CorrFile)
        
            fname = CorrFile(sp).name;
            indtmp1 = strfind(fname, keyword);
            spoolnum = sscanf(fname(indtmp1+length(keyword) : end), '%d', 1);
            indtmp2 = strfind(fname, cl(ii).lname);
            roinum = sscanf(fname(indtmp2+length(cl(ii).lname)+1 : end), '%d', 1);
            
            roiname = [keyword '_' num2str(spoolnum) '_' num2str(roinum)];
            
            CorrMatFile = load([pth '\' CorrFile(sp).name]);
            Correlation = double(CorrMatFile.CorrMat.Corr);
            density = double(CorrMatFile.CorrMat.density_1);
            Correlation(isnan(Correlation(:,2)), :) = [];
            Correlation(isinf(Correlation(:,2)), :) = [];

            if size(Correlation, 1) > 0
                para0 = [max(Correlation(:,2)), sqrt(2)*uncert(ii).xc(s), max(Correlation(:,2))/2, 100, 1];
                XDATA = Correlation(begin:end, 1);
                YDATA = Correlation(begin:end, 2);
                [parafit,~,Res,~,~,~,J] = lsqcurvefit(@autofit, para0, XDATA, YDATA, AutoLb, AutoUb, options_AC);
                YDATA_fit = autofit(parafit, XDATA);
                ChiSquare = sum((YDATA - YDATA_fit).^2./YDATA_fit);
                pval = 1 - chi2cdf(ChiSquare, size(YDATA, 1)-5);
                
                ci = nlparci(parafit, Res, 'jacobian', J);
                fprintf(fid_fitpara, fitparaFormat, roiname, parafit);
                fprintf(fid_fitparaSD, fitparaSDFormat, roiname, (ci(:,2) - ci(:,1))./2');
                
                if abs(parafit(2)) < abs(parafit(4))
                    Autofit(1) = 1/(parafit(1)*2*pi*parafit(2)^2);          % rho fluorophore 
                    Autofit(2) = parafit(2)/sqrt(2);                        % uncertainty  
                    Autofit(3) = sqrt(parafit(4)^2-parafit(2)^2)/sqrt(2);   % radii of cluster
                    Autofit(4) = parafit(3);                                % amplitude 
                    Autofit(5) = 2*Autofit(1)*Autofit(4)*pi*Autofit(3)^2;   % #proteins per cluster
                    Autofit(6) = Autofit(1)/Autofit(5);                     % cluster per nm^2
                else
                    Autofit(1) = 1/(parafit(3)*2*pi*parafit(4)^2);          % rho fluorophore 
                    Autofit(2) = parafit(4)/sqrt(2);                        % uncertainty  
                    Autofit(3) = sqrt(parafit(2)^2-parafit(4)^2)/sqrt(2);   % radii of cluster
                    Autofit(4) = parafit(1);                                % amplitude 
                    Autofit(5) = 2*Autofit(1)*Autofit(4)*pi*Autofit(3)^2;   % #proteins per cluster
                    Autofit(6) = Autofit(1)/Autofit(5);                     % cluster per nm^2
                end
                Autofit(7) = pval;
                fprintf(fid_fit, fitFormat, roiname, density, Autofit);
                
                h = 1;
                figure(h);
                semilogx(XDATA, YDATA, 'ko');
                hold on
                semilogx(XDATA, YDATA_fit,'r-','LineWidth',1.5);
                hold off
                saveas(h, [pth '\' fname(1:end-5) num2str(roinum) '.jpg']);
                clf(h);
                
                results = [Correlation(begin:end, 1:2), YDATA_fit];
                save([pth '\' fname(1:end-5) num2str(roinum) '_fit.txt'],'results','-ASCII');
            end
        end
        fclose(fid_fit);
        fclose(fid_fitpara);
        fclose(fid_fitparaSD);
    end
    
    CrossLb = [0, 0, 0, 0.95];
    CrossUb = [inf, inf, inf, inf];
    for ii = 1 : size(ColorMode, 2)
        for jj = ii+1 : size(ColorMode, 2)
            
            pth = [path0 '\' SampleList{s} '\PairCorrelation\' cl(ii).sname 'and' cl(jj).sname '\'];
            CrossFit = zeros(1, 3);
        
            fid_fit = fopen([path0 '\' SampleList{s} '\Cross' cl(ii).sname cl(jj).sname 'fit.txt'],'w');
            fprintf(fid_fit,'%-20s %-12s %-12s %-12s %-12s\r\n', 'ROI', 'Amp_unnorm', 'Amp_norm', 'Corrl', 'pval');
            fitFormat = '%-20s %-12.4E %-12.4E %-12.4E %-12.4E\r\n';
            
            fid_fitpara = fopen([path0 '\' SampleList{s} '\Cross' cl(ii).sname cl(jj).sname 'fitpara.txt'],'w');
            fprintf(fid_fitpara,'%-20s %-12s %-12s %-12s %-12s\r\n', 'ROI', 'Amp', 'xc', 'sigma', 'y0');
            fitparaFormat = '%-20s %-12.4E %-12.4E %-12.4E %-12.4E\r\n';
            
            fid_fitparaSD = fopen([path0 '\' SampleList{s} '\Cross' cl(ii).sname cl(jj).sname 'fitparaSD.txt'],'w');
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
                    
                    fprintf(fid_fit, fitFormat, roiname, density_1*density_2*CrossFit(1), CrossFit);
                    
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

function [calcu_ydata] = autofit(para, xdata)
%expression : y = y0 + a1*exp(-(x/r1)^2/2) + a2*exp(-(x/r2)^2/2)
% para = [a1 r1 a2 r2 y0] for auto fit
calcu_ydata = para(5) + para(1).*exp(-0.5.*xdata.^2./para(2).^2) +  para(3).*exp(-0.5*xdata.^2./para(4).^2);
end

function [calcu_ydata] = crossfit(para, xdata)
%expression : y = y0 + a*exp(-((x-c)/r)^2/2)+ a*exp(-((x+c)/r)^2/2)
% para = [a c r y0] for cross fit
calcu_ydata = para(4) + para(1).*exp(-0.5.*(xdata-para(2)).^2./para(3)^2) + para(1).*exp(-0.5.*(xdata+para(2)).^2./para(3)^2);
end