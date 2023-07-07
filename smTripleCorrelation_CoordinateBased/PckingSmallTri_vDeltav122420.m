clearvars
clc
fclose('all');
scrsz = get(groot, 'Screensize'); 
p = mfilename('fullpath');
[filepath, name, ~] = fileparts(p);
addpath(filepath);

kw = 'results';
% ======================== Setting up parameters ==========================
% ------------------------ perimeter threshold ----------------------------

thresh = 180; %nm

% ------------------------- Color Seq Info --------------------------------
TripleCalMode = 'BRG';
TriplelclMode = 'R';
switch strfind(TripleCalMode, TriplelclMode)
    case 1
        cross_col = 7;
        rho_col = 8;
    case 2
        cross_col = 5;
        rho_col = 9;
    case 3
        cross_col = 6;
        rho_col = 10;
end
        
%% =========================== Data Processing ============================
path0 = uigetdir('', 'Choose the Output Directory');
SampleList = dir(path0);
isSampleDir = [SampleList(:).isdir];
SampleNames = {SampleList(isSampleDir).name}';
SampleNames(ismember(SampleNames, {'.','..','map','Analysis','preview'})) = [];

for ii = 1 : numel(SampleNames)
    
    % ------------------------ Data Organizing ---------------------------
    TripleResultsFolder = ['TripleResults_' TripleCalMode '\Smooth_3_rng_100\'];
    files = dir([path0 '\' SampleNames{ii} '\' TripleResultsFolder '*.delta']);
    filenum = numel(files);
    filenames = {files(:).name}';
        
    fid_all = fopen([path0 '\' SampleNames{ii} '\' TripleResultsFolder 'PickedTris_delta_' num2str(thresh) '.txt'], 'w');
    fid_lcl = fopen([path0 '\' SampleNames{ii} '\' TripleResultsFolder 'PickedTris_delta_' num2str(thresh) '_' TriplelclMode '.preview'], 'w');
        
    for jj = 1 : filenum
            
        fname = filenames{jj}(length(kw)+1 : end-6);
            
        result_tmp = dlmread([path0 '\' SampleNames{ii} '\' TripleResultsFolder filenames{jj}], '', 1, 1);
        mask_peri = result_tmp(:, 1) + result_tmp(:, 2) + result_tmp(:, 3) > thresh; % perimeter thresholding
        %mask_neg = result_tmp(:, 4) < 0 | result_tmp(:, cross_col) <= 0; % negative thresholding
        mask_neg = result_tmp(:, 4) < 0 | result_tmp(:, 5) <= 0 | result_tmp(:, 6) <= 0 | result_tmp(:, 7) <= 0;
        mask = mask_peri | mask_neg;
        
        result_tmp(mask, :) = [];
        lcl_rho = result_tmp(:, 4) ./ result_tmp(:, cross_col) .* result_tmp(:, rho_col);
        
        save([path0 '\' SampleNames{ii} '\' TripleResultsFolder filenames{jj}(1:end-6) '_' num2str(thresh) '.flt'],'result_tmp','-ASCII');
        for kk = 1 : size(result_tmp, 1)
            fprintf(fid_all, '%-20s %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E %-12.4E\r\n', fname, result_tmp(kk, :));
            fprintf(fid_lcl, '%-20s %-12.4E %-12.4E %-12.4E\r\n', fname, result_tmp(kk, 4), result_tmp(kk, cross_col), lcl_rho(kk, 1));
        end
        
    end
    fclose(fid_all);
    fclose(fid_lcl);
        
    % ----------- preview the boxplots ----------------- 
    fid_lcl = fopen([path0 '\' SampleNames{ii} '\' TripleResultsFolder 'PickedTris_delta_' num2str(thresh) '_' TriplelclMode '.preview']);
    hh = textscan(fid_lcl, '%s %f %f %f');
    fclose(fid_lcl);
    
    gdata = log10(hh{4});
    cdata = hh{1};
    xdata = ones(size(cdata, 1), 1);
    for i = 1 : size(xdata, 1)-1
        if ~strcmp(cdata{i+1}, cdata{i})
            xdata(i+1 : end) = xdata(i+1 : end) + 1;
        end
    end
    h = figure('Position', scrsz);
    boxplot(gdata, cdata);
    hold on
    plot(xdata, gdata, 'ko');
    hold off
    saveas(h, [path0 '\' SampleNames{ii} '\' TripleResultsFolder 'PickedTris_delta_' num2str(thresh) '_' TriplelclMode '.jpg']);
    clf(h)
    close(h)
   
end

% ---------------------------- Outlier filtering --------------------------
numEntry = 0;
for ii = 1 : numel(SampleNames)
    
    TripleResultsFolder = ['TripleResults_' TripleCalMode '\Smooth_3_rng_100\'];
    fid = fopen([path0 '\' SampleNames{ii} '\' TripleResultsFolder 'PickedTris_delta_180_' TriplelclMode '.preview']);
    dummy = textscan(fid, '%s %f %f %f');
    fclose(fid);
    
    TF = isoutlier(log10(dummy{4}), 'quartiles');
    
    dummy{1}(TF) = [];
    dummy{2}(TF) = [];
    dummy{3}(TF) = [];
    dummy{4}(TF) = [];
    
    fid = fopen([path0 '\' SampleNames{ii} '\' TripleResultsFolder 'PickedTris_delta_180_' TriplelclMode '.qflted'], 'w');
    for jj = 1 : size(dummy{4}, 1)
        fprintf(fid, '%-20s %-12.4E %-12.4E %-12.4E\r\n', dummy{1}{jj}, dummy{2}(jj), dummy{3}(jj), dummy{4}(jj));
    end
    fclose(fid);
    numEntry = numEntry + size(dummy{4}, 1);

end

% -------------------------- overall boxplotting --------------------------
catagries{numEntry, 1} = [];
grps = zeros(numEntry, 1);
xdata = zeros(numEntry, 1);
flag = 0;
for ii = 1 : numel(SampleNames)
    
    TripleResultsFolder = ['TripleResults_' TripleCalMode '\Smooth_3_rng_100\'];
    fid = fopen([path0 '\' SampleNames{ii} '\' TripleResultsFolder 'PickedTris_delta_180_' TriplelclMode '.qflted']);
    dummy = textscan(fid, '%s %f %f %f');
    fclose(fid);
    
    catagries(flag+1 : flag+size(dummy{1}, 1)) = repmat(SampleNames(ii), size(dummy{1}, 1), 1);
    grps(flag+1 : flag+size(dummy{4}, 1)) = log10(dummy{4});
    xdata(flag+1 : flag+size(dummy{4}, 1)) = ii;
    
    flag = flag + size(dummy{4}, 1);

end
h = figure('Position', scrsz);
boxplot(grps, catagries);
xtickangle(45);
hold on
plot(xdata, grps, 'ko');
hold off
saveas(h, [path0 '\preview.jpg']);
clf(h);
close(h);