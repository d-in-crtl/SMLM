clearvars
clc
fclose('all');
addpath('C:\Users\Rothenberg\Documents\MATLAB\smTripleCorrelation\');
buffersz = 2000;
keyword_1 = 'spool_';
keyword_2 = 'triple_';

%% ==================== Perameter configurations ==========================
% ------------------------ seqential information --------------------------
%TripleCalMode  = ['BRG'; 'DBG'; 'DBR'; 'DRG'];
TripleCalMode = ['RBG'];
TripleCalMode  = cellstr(TripleCalMode);

%TripleDispMode = ['BGR'; 'BGD'; 'BRD'; 'RGD'];
TripleDispMode = ['BGR'];
TripleDispMode = cellstr(TripleDispMode);

% -------------------- triple correlation parameters ----------------------
rho_res = 5; % nm, usually triple correlation is calculated through 200*200 pixels with 5 nm/pixel
boxradii = [5, 5, 5]; % for smooth, in the unit of rho_res
dim = 500; % nm, correlation will be evaluated from 0 to /dim/ in all the 3 distances. 
dim = round(dim/rho_res);

% -------------------------- the norm mode --------------------------------
NormMode = 'const';
PlotMode = 'raw';
normDist = 200; % nm

% ----------------------- ploting graph parameters ------------------------
ColorCode = [229, 50, 56;  % ebay red
             134, 184, 23; % ebay green
             0, 100, 210;  % ebay blue
             245, 175, 2]; % ebay yellow
ColorCode = ColorCode./255;

BlackColor = [0 0 0];
BlackAlpha = 0.05;
ColorAlpha = 0.45;

magnification = [2000, 20, 50, 50];

%% =========================== Data reading ===============================
path0 = uigetdir('', 'Choose the Output Directory');
SampleList = dir(path0);
isSampleDir = [SampleList(:).isdir];
SampleNames = {SampleList(isSampleDir).name}';
SampleNames(ismember(SampleNames, {'.','..','map','Analysis','preview'})) = [];

%% =========================== Data processing ============================

for i = 1 : numel(SampleNames)
    for mm = 1 : numel(TripleCalMode)
        
        TriCalMode = TripleCalMode{mm}; 
        TriDispMode = TripleDispMode{mm}; 
    
        TripleResultsFolder = ['TripleFrcs_' TriCalMode '\'];
        
        % -------------- Identify color mode and corresponding cc ---------
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
        
        AutoFid_1 = fopen([path0 '\' SampleNames{i} '\Auto' ACdim_1 'fit.txt']);
        AutoFitResult_1 = textscan(AutoFid_1, '%s %f %f %f %f %f %f %f %f', 'HeaderLines', 1);
        AutoROI_1 = AutoFitResult_1{1};
        AutoRho_1 = AutoFitResult_1{2};
        fclose(AutoFid_1);
        
        files = dir([path0 '\' SampleNames{i} '\' TripleResultsFolder '*.mat']);
        filenames = {files(:).name}';
        filenum = numel(files) - sum(ismember(filenames, {'Avg_triple.mat','triple_trans.mat'}));
        filenames(ismember(filenames, {'Avg_triple.mat', 'triple_trans.mat'})) = [];
    
        results = zeros(buffersz, 5);
        Coordinates = zeros(buffersz, 6);
        Coordinates_norm = zeros(buffersz, 6);
        catamask = zeros(buffersz, 1); % recording its catagory
    
        NumofTris = 0;
    
        fid = fopen([path0 '\' SampleNames{i} '\' TripleResultsFolder 'results.localmax'],'w');
        fid_0 = fopen([path0 '\' SampleNames{i} '\' TripleResultsFolder 'results_cata_0.localmax'],'w');
        fid_1 = fopen([path0 '\' SampleNames{i} '\' TripleResultsFolder 'results_cata_1.localmax'],'w');
        fid_2 = fopen([path0 '\' SampleNames{i} '\' TripleResultsFolder 'results_cata_2.localmax'],'w');
        fid_3 = fopen([path0 '\' SampleNames{i} '\' TripleResultsFolder 'results_cata_3.localmax'],'w');
        fid_null = fopen([path0 '\' SampleNames{i} '\' TripleResultsFolder 'results_cata_null.localmax'],'w');
        
        fprintf(fid, '%-20s %-20s %-20s %-20s %-20s\r\n', ['DistOf' CCdim_1], ['DistOf' CCdim_2], ['DistOf' CCdim_3], ['FrcOf' TriCalMode(1)], ['RhoOf' TriCalMode(1)]);
        fprintf(fid_0, '%-20s %-20s %-20s %-20s %-20s\r\n', ['DistOf' CCdim_1], ['DistOf' CCdim_2], ['DistOf' CCdim_3], ['FrcOf' TriCalMode(1)], ['RhoOf' TriCalMode(1)]);
        fprintf(fid_1, '%-20s %-20s %-20s %-20s %-20s\r\n', ['DistOf' CCdim_1], ['DistOf' CCdim_2], ['DistOf' CCdim_3], ['FrcOf' TriCalMode(1)], ['RhoOf' TriCalMode(1)]);
        fprintf(fid_2, '%-20s %-20s %-20s %-20s %-20s\r\n', ['DistOf' CCdim_1], ['DistOf' CCdim_2], ['DistOf' CCdim_3], ['FrcOf' TriCalMode(1)], ['RhoOf' TriCalMode(1)]);
        fprintf(fid_3, '%-20s %-20s %-20s %-20s %-20s\r\n', ['DistOf' CCdim_1], ['DistOf' CCdim_2], ['DistOf' CCdim_3], ['FrcOf' TriCalMode(1)], ['RhoOf' TriCalMode(1)]);
        fprintf(fid_null, '%-20s %-20s %-20s %-20s %-20s\r\n', ['DistOf' CCdim_1], ['DistOf' CCdim_2], ['DistOf' CCdim_3], ['FrcOf' TriCalMode(1)], ['RhoOf' TriCalMode(1)]);
        
        triple_trans_avg = zeros(dim, dim, dim);
        for j = 1 : filenum
            
            fname = filenames{j};
            indtmp1 = strfind(fname, keyword_1);
            spoolnum = sscanf(fname(indtmp1+length(keyword_1) : end), '%d', 1);
            indtmp2 = strfind(fname, keyword_2);
            roinum = sscanf(fname(indtmp2+length(keyword_2) : end), '%d', 1);
            
            fname = filenames{j}(1:end-4);
            triple = load([path0 '\' SampleNames{i} '\' TripleResultsFolder  filenames{j}]);
            fieldname = ['TripleResult_' TriCalMode];
            triplemat = triple.(fieldname);
            
            idx_AC_1 = find(strcmp(AutoROI_1, ['spool__' num2str(spoolnum) '_' num2str(roinum)]));
            Rho_AC_1 = AutoRho_1(idx_AC_1);
            
            % ----------------- processing raw triple cube --------------------
            
            [triple_trans, results_tmp, coordinates_tmp, coordinates_tmp_norm] = TripleRead_Frc(triplemat, boxradii, dim, rho_res, NormMode, normDist, TriCalMode, TriDispMode);
            triple_trans_avg = triple_trans_avg + triple_trans;
            
            % ------------------------ Catagorize -----------------------------
            if~isempty(coordinates_tmp)
                
                x1 = coordinates_tmp(:,1);
                y1 = coordinates_tmp(:,2);
                x2 = coordinates_tmp(:,3);
                y2 = coordinates_tmp(:,4);
                x3 = coordinates_tmp(:,5);
                y3 = coordinates_tmp(:,6);
                
                cata_1 = (x3 > -1.5.*(x2-x1) & x3 < -0.5.*(x2-x1)) & (abs(y3) <= (x2-x1)./2 | abs(y3) <= (abs(x3)-0.5.*(x2-x1)).*tan(pi/4));
                cata_2 = (x3 > -0.5.*(x2-x1) & x3 <  0.5.*(x2-x1)) & (abs(y3) <= (x2-x1)./2);
                cata_3 = (x3 >  0.5.*(x2-x1) & x3 <  1.5.*(x2-x1)) & (abs(y3) <= (x2-x1)./2 | abs(y3) <= (abs(x3)-0.5.*(x2-x1)).*tan(pi/4));
                cata_0 = (~cata_1 & ~cata_2 & ~cata_3) & (results_tmp(:, 1) < 100 & results_tmp(:, 2) < 100 & results_tmp(:, 3) < 100);
                
                NumofResults = size(results_tmp, 1);
                mask_tmp = -1 .* ones(NumofResults, 1);
                mask_tmp(cata_0, :) = 0;
                mask_tmp(cata_1, :) = 1;
                mask_tmp(cata_2, :) = 2;
                mask_tmp(cata_3, :) = 3;
        
                % -------------------- Plotting and Preview -----------------------
                h = TriplePlot_sn(triple_trans, dim, rho_res, results_tmp, coordinates_tmp, mask_tmp, magnification(mm), TriCalMode, TriDispMode, ColorCode, BlackColor, BlackAlpha, ColorAlpha);
                saveas(h,[path0 '\' SampleNames{i} '\' TripleResultsFolder  fname '.jpg'],'jpg');
                clf(h); 
                close(h);
        
                % --------------------- Collecting results ------------------------
                results(NumofTris+1 : NumofTris+NumofResults, 1:4) = results_tmp;
                results(NumofTris+1 : NumofTris+NumofResults, 5) = Rho_AC_1 .* ones(NumofResults, 1);
                
                Coordinates(NumofTris+1 : NumofTris+NumofResults, :) = coordinates_tmp;
                Coordinates_norm(NumofTris+1 : NumofTris+NumofResults, :) = coordinates_tmp_norm;
                catamask(NumofTris+1 : NumofTris+NumofResults, :) = mask_tmp;
        
                NumofTris = NumofTris + NumofResults;
            end
        end
        triple_trans_avg = triple_trans_avg / filenum;
        save([path0 '\' SampleNames{i} '\' TripleResultsFolder 'triple_trans.mat'],'triple_trans_avg');
        
        results((NumofTris+1):end, :) = [];
        Coordinates((NumofTris+1):end, :) = [];
        Coordinates_norm((NumofTris+1):end, :) = [];
        catamask((NumofTris+1):end, :) = [];
    
        % ---------------- Plotting and printing all the results --------------
        [h_allbut0, h0, h1, h2, h3] = TriplePlot_All(results, Coordinates, Coordinates_norm, catamask, magnification(mm), TriDispMode, ColorCode, BlackColor, BlackAlpha, ColorAlpha, PlotMode);
        clear Coordinates_norm Coordinates;
        
        print(h_allbut0,[path0 '\' SampleNames{i} '\' TripleResultsFolder 'stat_total.pdf'],'-dpdf','-bestfit', '-painters');
        clf(h_allbut0); close(h_allbut0);
        print(h0,[path0 '\' SampleNames{i} '\' TripleResultsFolder 'stat_overlap.pdf'],'-dpdf','-bestfit', '-painters');
        clf(h0); close(h0);
        print(h1,[path0 '\' SampleNames{i} '\' TripleResultsFolder 'stat_cata1.pdf'],'-dpdf','-bestfit', '-painters');
        clf(h1); close(h1);
        print(h2,[path0 '\' SampleNames{i} '\' TripleResultsFolder 'stat_cata2.pdf'],'-dpdf','-bestfit', '-painters');
        clf(h2); close(h2);
        print(h3,[path0 '\' SampleNames{i} '\' TripleResultsFolder 'stat_cata3.pdf'],'-dpdf','-bestfit', '-painters');
        clf(h3); close(h3);
    
        % ------------------- Saving results tables ---------------------------
        results_cata_null = results(catamask(:, 1) == -1, :);
        results_cata_0 = results(catamask(:, 1) == 0, :);
        results_cata_1 = results(catamask(:, 1) == 1, :);
        results_cata_2 = results(catamask(:, 1) == 2, :);
        results_cata_3 = results(catamask(:, 1) == 3, :);
        clear catamask
    
        for j = 1 : size(results, 1)
            fprintf(fid, '%-12.4E %-12.4E %-12.4E %-12.4E %-12.4E\r\n', results(j, :));
        end
        clear results
        fclose(fid);
    
        for j = 1 : size(results_cata_0, 1)
            fprintf(fid_0, '%-12.4E %-12.4E %-12.4E %-12.4E %-12.4E\r\n', results_cata_0(j, :));
        end
        clear results_cata_0
        fclose(fid_0);
    
        for j = 1 : size(results_cata_1, 1)
            fprintf(fid_1, '%-12.4E %-12.4E %-12.4E %-12.4E %-12.4E\r\n', results_cata_1(j, :));
        end
        clear results_cata_1
        fclose(fid_1);
    
        for j = 1 : size(results_cata_2, 1)
            fprintf(fid_2, '%-12.4E %-12.4E %-12.4E %-12.4E %-12.4E\r\n', results_cata_2(j, :));
        end
        clear results_cata_2
        fclose(fid_2);
    
        for j = 1 : size(results_cata_3, 1)
            fprintf(fid_3, '%-12.4E %-12.4E %-12.4E %-12.4E %-12.4E\r\n', results_cata_3(j, :));
        end
        clear results_cata_3
        fclose(fid_3);
        
        for j = 1 : size(results_cata_null, 1)
            fprintf(fid_null, '%-12.4E %-12.4E %-12.4E %-12.4E %-12.4E\r\n', results_cata_null(j, :));
        end
        clear results_cata_null
        fclose(fid_null);
    
    end
end