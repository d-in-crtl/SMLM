clearvars
clc
fclose('all');

p = mfilename('fullpath');
[filepath, name, ~] = fileparts(p);
addpath(filepath);

% ========================== Color info ===================================
Color = 'RED'; % Only Supports RED, GREEN, BLUE, and DarkRed
% ====================== DBSCAN parameters ================================
eps = 15; % nm
min_points = 1;

% =================== Camera and Reconstruction Info ======================
ReconMode = 'MFA_new';
Camera_pxsz = 65.39; % nm
Image_pxsz = 65.39; % nm. This is the pixel size of the image from which you draw the ROI. 
                 
                 % For instance, 
                 % if the ROI was obtained from .frm1, then Image_pxsz = 73.3 nm;
                 % if the ROI was obtained from Reonstruction.tif, then the Image_pxsz = 10 nm.
                 
                 % ROI can be drawn freehand, but large size of ROI might cause failue in computer RAM.
% =========================================================================

path0 = uigetdir('', 'Choose the Output Directory');
SampleList = dir(path0);
isSampleDir = [SampleList(:).isdir];
SampleNames = {SampleList(isSampleDir).name}';
SampleNames(ismember(SampleNames, {'.','..','map','Analysis','preview'})) = [];

for s = 1 : numel(SampleNames)   
    TableResults = dir([path0 '\' SampleNames{s} '\*_' Color '.result']);
    filesz = [TableResults(:).bytes];
    TableNames = {TableResults(filesz > 1000).name}';
    
    for sp = 1 : numel(TableNames)   
        FileName = TableNames{sp};
        FullFileName = [path0 '\' SampleNames{s} '\' FileName];
        report_table = dlmread(FullFileName);
        
        ZipName = [FileName(1:end-length(Color)-length('.result')) 'roi.zip'];
        FullZipName = [path0 '\' SampleNames{s} '\' ZipName];
        
        [roiObj] = readroi(FullZipName);
        roi_coors = getROIcoors(report_table, roiObj, Camera_pxsz, Image_pxsz, ReconMode); % output coordinates' unit is 'nm'
        
        for r = 1 : numel(roiObj)    
            fid_Clusters = fopen([path0 '\' SampleNames{s} '\' FileName(1:end-length('.result')) '_roi_' num2str(r) '.Cluster'], 'w');
            fprintf(fid_Clusters, '%s %s %s %s %s %s\r\n', 'ClusterID', 'NumSpots', 'Area_Ellipse', 'density', 'center_x', 'center_y');
            
            fid_Coors = fopen([path0 '\' SampleNames{s} '\' FileName(1:end-length('.result')) '_roi_' num2str(r) '.Coors'], 'w');
            fprintf(fid_Coors, '%s %s %s\r\n', 'ClusterID', 'x (nm)', 'y (nm)');
            
            disp([num2str(s) ': ' SampleNames{s} ' -' FileName(1:end-length('.result')) '-roi-' num2str(r)]);
            
            tic; 
            
            x = roi_coors(r).x;
            y = roi_coors(r).y;
            points = [x, y];
            
            if isempty(points)
                continue
            end
            [coors_table, Cluster_List, Cluster_Figure] = mst_dbscan(points, eps, min_points);
            
            
            fprintf(fid_Clusters, '%f %f %f %f %f %f\r\n', Cluster_List');
            fprintf(fid_Coors, '%f %f %f\r\n', coors_table');
            savefig(Cluster_Figure, [path0 '\' SampleNames{s} '\' FileName(1:end-length('.result')) '_roi_' num2str(r) '.fig']);
            saveas(Cluster_Figure, [path0 '\' SampleNames{s} '\' FileName(1:end-length('.result')) '_roi_' num2str(r) '.jpg']);
            close(Cluster_Figure);
            fclose(fid_Clusters);
            fclose(fid_Coors);
            
            toc;
        end
    end
end
