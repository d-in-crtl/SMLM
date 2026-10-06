clearvars
clc
fclose('all');

p = mfilename('fullpath');
[filepath, name, ~] = fileparts(p);
addpath(filepath);

% ========================== Color info ===================================
Color = 'RED'; % Only Supports RED, GREEN, BLUE, and DarkRed


bin_step = 1;
bin_max = 500;
bin_sz = (5 : bin_step : bin_max); % bin size of number of molecules per cluster
% =========================================================================

path0 = uigetdir('', 'Choose the Output Directory');
SampleList = dir(path0);
isSampleDir = [SampleList(:).isdir];
SampleNames = {SampleList(isSampleDir).name}';
SampleNames(ismember(SampleNames, {'.','..','map','Analysis','preview'})) = [];

sz_hist = zeros(2, size(bin_sz, 2) - 1);
sz_hist(1, :) = bin_sz(2 : end) - bin_step / 2;
for s = 1 : numel(SampleNames)
    
    ClusterFileList = dir([path0 '\' SampleNames{s} '\*_' Color '_*.Cluster']);
    ClusterFileNames = {ClusterFileList(:).name}';
    
    for i = 1 : numel(ClusterFileNames)
        
        FileName = ClusterFileNames{i};
        fid = fopen([path0 '\' SampleNames{s} '\' FileName]);
        clusters = textscan(fid, '%f %f %f %f %f %f', 'HeaderLines', 1);
        fclose(fid);
        
        sz = clusters{2};
        sz_hist(2, :) = sz_hist(2, :) + histcounts(sz, bin_sz);
    
    end
   
end
fid = fopen([path0 '\szHistOf' Color '.hist'], 'w');
fprintf(fid, '%f %f\r\n', sz_hist);
fclose(fid);
