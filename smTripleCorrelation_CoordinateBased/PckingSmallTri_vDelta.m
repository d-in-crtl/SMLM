clearvars
clc
fclose('all');
addpath('C:\Users\Rothenberg\Documents\MATLAB\smTripleCorrelation\');
buffersz = 2000;
thresh = 180; %nm
% ======================== Setting up parameters ==========================
% ------------------------- Color Seq Info --------------------------------
%TripleCalMode  = ['BRG'; 'DBG'; 'DBR'; 'DRG'];
TripleCalMode = ['BRG'];
TripleCalMode  = cellstr(TripleCalMode);

%TripleDistMode = ['BG'; 'BG'; 'BR'; 'RG'];
TripleDistMode = ['BG'];
TripleDistMode = cellstr(TripleDistMode);
abc = [1, 2, 3];

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
        TriDistMode = TripleDistMode{mm}; 
        seq = zeros(1,2);
        seq(1) = find(TriCalMode == TriDistMode(1));
        seq(2) = find(TriCalMode == TriDistMode(2));
        sequence = min(seq(1), seq(2))*10 + max(seq(1), seq(2));
        switch sequence
            case 13
                seqmode = 1;
            case 12
                seqmode = 2;
            case 23
                seqmode = 3;
        end
        % ------------------------ Data Reading ---------------------------
        TripleResultsFolder = ['TripleResults_' TriCalMode '\Smooth_3_rng_100\'];
        files = dir([path0 '\' SampleNames{ii} '\' TripleResultsFolder '*.delta']);
        filenum = numel(files);
        filenames = {files(:).name}';
        
        result_localmax = zeros(buffersz, 10);
        NumofTri = 0;
        for jj = 1 : filenum
            result_tmp = dlmread([path0 '\' SampleNames{ii} '\' TripleResultsFolder filenames{jj}], '', 1, 1);
            mask_tmp = result_tmp(:, 1) + result_tmp(:, 2) + result_tmp(:, 3) > thresh;
            result_tmp(mask_tmp, :) = [];
            save([path0 '\' SampleNames{ii} '\' TripleResultsFolder filenames{jj}(1:end-6) '_' num2str(thresh) '.flt'],'result_tmp','-ASCII');
            
            result_localmax(NumofTri+1 : NumofTri+size(result_tmp, 1), :) = result_tmp;
            NumofTri = NumofTri + size(result_tmp, 1);
        end
        result_localmax(NumofTri+1 : end, :) = [];
        save([path0 '\' SampleNames{ii} '\' TripleResultsFolder 'PickedTris_delta_' num2str(thresh) '.txt'],'result_localmax','-ASCII');
        
        %print(h_allbut0,[path0 '\' SampleNames{ii} '\' TripleResultsFolder 'stat_total.pdf'],'-dpdf','-bestfit', '-painters');
        %clf(h_allbut0); close(h_allbut0);
        
        
    end
end