clearvars
clc
fclose('all');
addpath('C:\Users\Rothenberg\Documents\MATLAB\smTripleCorrelation\');
buffersz = 2000;

% ======================== Setting up parameters ==========================
% ------------------------- Color Seq Info --------------------------------
TripleCalMode  = ['BRG'; 'DBG'; 'DBR'; 'DRG'];
%TripleCalMode = ['RGB'];
TripleCalMode  = cellstr(TripleCalMode);

TripleDispMode = ['RGB'; 'BGD'; 'RBD'; 'RGD'];
%TripleDispMode = ['RGB'];
TripleDispMode = cellstr(TripleDispMode);

% -------------------------- the norm mode --------------------------------
NormMode = 'const';
PlotMode = 'norm';
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

magnification = [150, 35, 50, 50];

%% =========================== Data Processing ============================
path0 = uigetdir('', 'Choose the Output Directory');
SampleList = dir(path0);
isSampleDir = [SampleList(:).isdir];
SampleNames = {SampleList(isSampleDir).name}';
SampleNames(ismember(SampleNames, {'.','..','map','Analysis','preview'})) = [];

for ii = 1 : numel(SampleNames)
    for mm = 1 : numel(TripleCalMode)
        
        % ------------------------ DispSeq Info ---------------------------
        TriCalMode = TripleCalMode{mm}; 
        TriDispMode = TripleDispMode{mm}; 
        seq = zeros(1,3);
        seq(1) = find(TriCalMode == TriDispMode(1));
        seq(2) = find(TriCalMode == TriDispMode(2));
        seq(3) = find(TriCalMode == TriDispMode(3));
        seqmode = seq(1)*100 + seq(2)*10 + seq(3);
        
        % ------------------------ Data Reading ---------------------------
        TripleResultsFolder = ['TripleResults_' TriCalMode '\'];
        files = dir([path0 '\' SampleNames{ii} '\' TripleResultsFolder '*.localmax']);
        filenum = numel(files);
        filenames = {files(:).name}';
        
        result_localmax = zeros(buffersz, 4);
        NumofTri = 0;
        for jj = 1 : filenum
            result_tmp = dlmread([path0 '\' SampleNames{ii} '\' TripleResultsFolder filenames{jj}], '', 1, 0);
            result_localmax(NumofTri+1 : NumofTri+size(result_tmp, 1), :) = result_tmp;
            NumofTri = NumofTri + size(result_tmp, 1);
        end
        result_localmax(NumofTri+1 : end, :) = [];
        
        indx_12 = result_localmax(:, 2);
        indx_13 = result_localmax(:, 1);
        indx_23 = result_localmax(:, 3);
        
        % ----------------- Translating Edges to Coordinates --------------
        [coordinates, coordinates_norm] = Edge2Coor(indx_13, indx_12, indx_23, seqmode, NormMode, normDist);
        
        % ------------------------- Catagorize ----------------------------
        if~isempty(coordinates)
                
            x1 = coordinates(:, 1);
            y1 = coordinates(:, 2);
            x2 = coordinates(:, 3);
            y2 = coordinates(:, 4);
            x3 = coordinates(:, 5);
            y3 = coordinates(:, 6);
                
            cata_1 = (x3 > -1.5.*(x2-x1) & x3 < -0.5.*(x2-x1)) & (abs(y3) <= (x2-x1)./2 | abs(y3) <= (abs(x3)-0.5.*(x2-x1)).*tan(pi/4));
            cata_2 = (x3 > -0.5.*(x2-x1) & x3 <  0.5.*(x2-x1)) & (abs(y3) <= (x2-x1)./2);
            cata_3 = (x3 >  0.5.*(x2-x1) & x3 <  1.5.*(x2-x1)) & (abs(y3) <= (x2-x1)./2 | abs(y3) <= (abs(x3)-0.5.*(x2-x1)).*tan(pi/4));
            cata_0 = (~cata_1 & ~cata_2 & ~cata_3) & (result_localmax(:, 1) < 80 & result_localmax(:, 2) < 80 & result_localmax(:, 3) < 80);
                
            NumofResults = size(coordinates, 1);
            catamask = -1 .* ones(NumofResults, 1);
            catamask(cata_0, :) = 0;
            catamask(cata_1, :) = 1;
            catamask(cata_2, :) = 2;
            catamask(cata_3, :) = 3;
        end
        
        % ----------------------- Plotting and Printing -------------------
        [h_allbut0, h0, h1, h2, h3] = TriplePlot_All(result_localmax, coordinates, coordinates_norm, catamask, magnification(mm), TriDispMode, ColorCode, BlackColor, BlackAlpha, ColorAlpha, PlotMode);
        
        print(h_allbut0,[path0 '\' SampleNames{ii} '\' TripleResultsFolder 'stat_total.pdf'],'-dpdf','-bestfit', '-painters');
        clf(h_allbut0); close(h_allbut0);
        print(h0,[path0 '\' SampleNames{ii} '\' TripleResultsFolder 'stat_overlap.pdf'],'-dpdf','-bestfit', '-painters');
        clf(h0); close(h0);
        print(h1,[path0 '\' SampleNames{ii} '\' TripleResultsFolder 'stat_cata1.pdf'],'-dpdf','-bestfit', '-painters');
        clf(h1); close(h1);
        print(h2,[path0 '\' SampleNames{ii} '\' TripleResultsFolder 'stat_cata2.pdf'],'-dpdf','-bestfit', '-painters');
        clf(h2); close(h2);
        print(h3,[path0 '\' SampleNames{ii} '\' TripleResultsFolder 'stat_cata3.pdf'],'-dpdf','-bestfit', '-painters');
        clf(h3); close(h3);
        
    end
end