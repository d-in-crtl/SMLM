clearvars
clc
fclose('all');

p = mfilename('fullpath');
[filepath, name, ~] = fileparts(p);
addpath(filepath); 

thresh = 400; % nm, distance threshold
%The threshold is to define whether AB is associated. 
%For example, from spool_i_colorA_roi.NND, you get the NND for EACH colorA cluster. 
%If this NND is larger than the threshold, then this colorA cluster does not form AB, vice versa otherwise.

Path0 = 'Z:\homes\guptad07\DBScan_ESCO2\TimelessAtG4andForks\070522\TimelessatForks\Trial_15g5\Trial_15r20';
SampleLists = dir(Path0);
isSampleDir = [SampleLists(:).isdir];
SampleNames = {SampleLists(isSampleDir).name}';
SampleNames(ismember(SampleNames, {'.','..','map','Analysis','preview'})) = [];

kw_1 = 'spool_';
kw_2 = '_roi_';

%% CF 2021-07-17
AllDist = [];
AllDistRndm = [];
%%

for s = 1 : numel(SampleNames)
    
    Path1 = [Path0 '\' SampleNames{s} '\NND_REDvsGREEN'];
    NNDLists = dir([Path1 '\spool*_GREEN_*.NND']);
    NNDNames = {NNDLists(:).name}';
    %keyboard
    fid_sv = fopen([Path1 '\NNDFracGREEN_thresh_' num2str(thresh) '.txt'], 'w'); % NNDfrac_colorA is the fraction of colorA that form AB, namely AB/A.  
    fprintf(fid_sv, '%s %s %s %s\r\n', 'ROI', 'num_GREEN', 'Frac_NND2RED', 'Frac_NND2RED_RND');
    
    for n = 1 : numel(NNDNames)
        
        fname = [Path1 '\' NNDNames{n}];
        
        indtmp1 = strfind(fname, kw_1);
        sp_num = sscanf(fname(indtmp1+length(kw_1) : end), '%d', 1);
        
        indtmp2 = strfind(fname, kw_2);
        roi_num = sscanf(fname(indtmp2+length(kw_2) : end), '%d', 1);
        
        ROI = [kw_1 num2str(sp_num) kw_2 num2str(roi_num)];
        
        fid_rd = fopen([Path1 '\' NNDNames{n}]);
        dummy = textscan(fid_rd, '%f %f %f %f', 'HeaderLines', 1); % CF changed from dummy = textscan(fid_rd, '%f %f %f', 'HeaderLines', 1);
        fclose(fid_rd);
%         keyboard
        dist = dummy{3}; % CF changed from dist = dummy{2}
        dist_rnd = dummy{4}; % CF changed from dist_rnd = dummy{3};
%         keyboard
        frac = sum(dist <= thresh) / size(dist, 1);
        frac_rnd = sum(dist_rnd <= thresh) / size(dist, 1);
        fprintf(fid_sv, '%s %d %f %f\r\n', ROI, size(dist, 1), frac, frac_rnd);
        AllDist = [AllDist; dist];
        AllDistRndm = [AllDistRndm; dist_rnd];
    end
    fclose(fid_sv);
end
