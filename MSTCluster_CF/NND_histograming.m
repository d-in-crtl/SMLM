clearvars
clc
fclose('all');

p = mfilename('fullpath');
[filepath, name, ~] = fileparts(p);
addpath(filepath);

NND_binsz = 1; % nm
NND_maxrange = 1000; % nm

NND_bins = (0 : NND_binsz : NND_maxrange);

Path0 = 'Z:\homes\guptad07\DBScan_ESCO2\TimelessAtG4andForks\061422\TimelessAtForks_070223\Trail_15r5\Trail_15g20';
SampleLists = dir(Path0);
isSampleDir = [SampleLists(:).isdir];
SampleNames = {SampleLists(isSampleDir).name}';
SampleNames(ismember(SampleNames, {'.','..','map','Analysis','preview'})) = [];

for s = 1 : numel(SampleNames)
    
    NND_hist = zeros(1, length(NND_bins) - 1);
    NND_rnd_hist = zeros(1, length(NND_bins) - 1);
    
    Path1 = [Path0 '\' SampleNames{s} '\NND_REDvsGREEN'];
    NNDLists = dir([Path1 '\spool*RED*.NND']);
    NNDNames = {NNDLists(:).name}';
    
    fid_sv = fopen([Path1 '\NND_hists_GREEN2RED.hist'], 'w'); % NND from red to each dark red cluster
    fprintf(fid_sv, '%12s %12s %12s\r\n', 'NND', 'hist', 'hist_rnd');
    
    for n = 1 : numel(NNDNames)
        
        fid_rd = fopen([Path1 '\' NNDNames{n}]);
        dummy = textscan(fid_rd, '%f %f %f', 'HeaderLines', 1);
        fclose(fid_rd);
        dist = dummy{2};
        dist_rnd = dummy{3};
        
        NND_hist = NND_hist + histcounts(dist, NND_bins);
        NND_rnd_hist = NND_rnd_hist + histcounts(dist_rnd, NND_bins);
        
    end
    
    result = [NND_bins(2:end) - NND_binsz/2; NND_hist; NND_rnd_hist];
    fprintf(fid_sv, '%f %f %f\r\n', result);
    fclose(fid_sv);
end
