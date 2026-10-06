clearvars
clc
fclose('all');

p = mfilename('fullpath');
[filepath, name, ~] = fileparts(p);
addpath(filepath);

% ========================== Color info ===================================
Colors = ['GREEN';
          'RED  '];% Only Supports RED, GREEN, BLUE, and DarkRed
Color = cellstr(Colors);

% =================== Camera and Reconstruction Info ======================
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
    SpoolZips = dir([path0 '\' SampleNames{s} '\*_roi.zip']);
    
    for sp = 1 : numel(SpoolZips)    
        ZipName = SpoolZips(sp).name;
        SpoolName = ZipName(1:end-8);
        fname1 = [path0 '\' SampleNames{s} '\' SpoolName '_' Color{1} '.result'];
        fname2 = [path0 '\' SampleNames{s} '\' SpoolName '_' Color{2} '.result'];
        if ~isfile(fname1) || ~isfile(fname2)
            continue
        end
        FullZipName = [path0 '\' SampleNames{s} '\' ZipName];
        [roiObj] = readroi(FullZipName);
        
        for r = 1 : numel(roiObj)    
            range = pdist2([roiObj(r).bdry(:, 1).*Image_pxsz, roiObj(r).bdry(:, 2).*Image_pxsz], [roiObj(r).bdry(:, 1).*Image_pxsz, roiObj(r).bdry(:, 2).*Image_pxsz]);
            MaxDist = max(range(:));
            
            % reading the coorfile that classified into different cluster ID
            Coor_fid_1 = fopen([path0 '\' SampleNames{s} '\spool_' num2str(sp) '_' Color{1} '_roi_' num2str(r) '.Coors']);
            CoorFullFile_1 = textscan(Coor_fid_1, '%f %f %f', 'HeaderLines', 1);
            fclose(Coor_fid_1);
            
            Cluster_fid_1 = fopen([path0 '\' SampleNames{s} '\spool_' num2str(sp) '_' Color{1} '_roi_' num2str(r) '.Cluster']);
            ClusterFullFile_1 = textscan(Cluster_fid_1, '%f %f %f %f %f %f', 'HeaderLines', 1);
            fclose(Cluster_fid_1);
            ClusterID_1 = ClusterFullFile_1{1};
%             keyboard
            Coor_fid_2 = fopen([path0 '\' SampleNames{s} '\spool_' num2str(sp) '_' Color{2} '_roi_' num2str(r) '.Coors']);
            CoorFullFile_2 = textscan(Coor_fid_2, '%f %f %f', 'HeaderLines', 1);
            fclose(Coor_fid_2);
            
            Cluster_fid_2 = fopen([path0 '\' SampleNames{s} '\spool_' num2str(sp) '_' Color{2} '_roi_' num2str(r) '.Cluster']);
            ClusterFullFile_2 = textscan(Cluster_fid_2, '%f %f %f %f %f %f', 'HeaderLines', 1);
            fclose(Cluster_fid_2);
            ClusterID_2 = ClusterFullFile_2{1};
            
            if isempty(ClusterID_1) || isempty(ClusterID_2)
                continue
            end
            
            CoorFile_1 = [CoorFullFile_1{1}, CoorFullFile_1{2}, CoorFullFile_1{3}];
            CoorFile_2 = [CoorFullFile_2{1}, CoorFullFile_2{2}, CoorFullFile_2{3}];
            
            % number of clusters for both channels
            numCluster_1 = size(ClusterID_1, 1);
            numCluster_2 = size(ClusterID_2, 1);
            
            disp([num2str(s) ': ' SampleNames{s} '-' SpoolName '-roi-' num2str(r)]);
            disp([Color{1} '-' num2str(numCluster_1) ', ' Color{2} '-' num2str(numCluster_2)]);
            
            tic; 
            
            % Computing NND and NND after randomization
            if ~isfolder([path0 '\' SampleNames{s} '\NND_' Color{1} 'vs' Color{2}])
                mkdir([path0 '\' SampleNames{s} '\NND_' Color{1} 'vs' Color{2}]);
            end
            
            fid_NND_1 = fopen([path0 '\' SampleNames{s} '\NND_' Color{1} 'vs' Color{2} '\' SpoolName '_' Color{1} '_roi_' num2str(r) '.NND'], 'w');
            fprintf(fid_NND_1, '%s %s %s %s\r\n', 'ClusterID', 'Contents', ['NND2' Color{2}], ['NND2' Color{2} '_RND']);
            dumidx = ismember(uint32(CoorFile_2(:, 1)), uint32(ClusterID_2));
            for ii = 1 : numCluster_1
                % cluster coordinates
                idx = uint32(CoorFile_1(:, 1)) == ClusterID_1(ii);
                numMol = sum(idx);
                tmp_x_ini = CoorFile_1(idx, 2);
                tmp_y_ini = CoorFile_1(idx, 3);
                center_x_ini = mean(tmp_x_ini);
                center_y_ini = mean(tmp_y_ini);
                
                % randomazing cluster position
                tmp_r = MaxDist * rand;
                tmp_theta = 2*pi * rand;
                center_x_rnd = center_x_ini + tmp_r * cos(tmp_theta);
                center_y_rnd = center_y_ini + tmp_r * sin(tmp_theta);
                while isempty(inpolygon(center_x_rnd, center_y_rnd, roiObj(r).bdry(:, 1).*Image_pxsz, roiObj(r).bdry(:, 2).*Image_pxsz))
                    tmp_r = MaxDist * rand;
                    tmp_theta = 2*pi * rand;
                    center_x_rnd = mean(tmp_x_ini) + tmp_r * cos(tmp_theta);
                    center_y_rnd = mean(tmp_y_ini) + tmp_r * sin(tmp_theta);
                end
                % randomazing cluster orientation
                tmp_rot = 2*pi*rand;
                tmp_x_rnd = (tmp_x_ini + tmp_r * cos(tmp_theta) - center_x_rnd) .* cos(tmp_rot) - (tmp_y_ini + tmp_r * sin(tmp_theta) - center_y_rnd) .* sin(tmp_rot) + center_x_rnd;
                tmp_y_rnd = (tmp_x_ini + tmp_r * cos(tmp_theta) - center_x_rnd) .* sin(tmp_rot) + (tmp_y_ini + tmp_r * sin(tmp_theta) - center_y_rnd) .* cos(tmp_rot) + center_y_rnd;
                
                % computing NND
                D2_ini = pdist2([tmp_x_ini, tmp_y_ini], [CoorFile_2(dumidx, 2), CoorFile_2(dumidx, 3)]);
                NND2_ini = min(D2_ini(:));
                
                % computing randomized NND
                D2_rnd = pdist2([tmp_x_rnd, tmp_y_rnd], [CoorFile_2(dumidx, 2), CoorFile_2(dumidx, 3)]);
                NND2_rnd = min(D2_rnd(:));
                
                fprintf(fid_NND_1, '%f %f %f %f\r\n', ii, numMol, NND2_ini, NND2_rnd);
            end
            fclose(fid_NND_1);
            
            fid_NND_2 = fopen([path0 '\' SampleNames{s} '\NND_' Color{1} 'vs' Color{2} '\' SpoolName '_' Color{2} '_roi_' num2str(r) '.NND'], 'w');
            fprintf(fid_NND_2, '%s %s %s\r\n', 'ClusterID', ['NND2' Color{1}], ['NND2' Color{1} '_RND']);
            dumidx = ismember(uint32(CoorFile_1(:, 1)), uint32(ClusterID_1));
            for ii = 1 : numCluster_2
                % cluster coordinates
                idx = uint32(CoorFile_2(:, 1)) == ClusterID_2(ii);
                numMol = sum(idx);
                tmp_x_ini = CoorFile_2(idx, 2);
                tmp_y_ini = CoorFile_2(idx, 3);
                center_x_ini = mean(tmp_x_ini);
                center_y_ini = mean(tmp_y_ini);
                
                % randomazing cluster position
                tmp_r = MaxDist * rand;
                tmp_theta = 2*pi * rand;
                center_x_rnd = center_x_ini + tmp_r * cos(tmp_theta);
                center_y_rnd = center_y_ini + tmp_r * sin(tmp_theta);
                while isempty(inpolygon(center_x_rnd, center_y_rnd, roiObj(r).bdry(:, 1).*Image_pxsz, roiObj(r).bdry(:, 2).*Image_pxsz))
                    tmp_r = MaxDist * rand;
                    tmp_theta = 2*pi * rand;
                    center_x_rnd = mean(tmp_x_ini) + tmp_r * cos(tmp_theta);
                    center_y_rnd = mean(tmp_y_ini) + tmp_r * sin(tmp_theta);
                end
                % randomazing cluster orientation
                tmp_rot = 2*pi*rand;
                tmp_x_rnd = (tmp_x_ini + tmp_r * cos(tmp_theta) - center_x_rnd) .* cos(tmp_rot) - (tmp_y_ini + tmp_r * sin(tmp_theta) - center_y_rnd) .* sin(tmp_rot) + center_x_rnd;
                tmp_y_rnd = (tmp_x_ini + tmp_r * cos(tmp_theta) - center_x_rnd) .* sin(tmp_rot) + (tmp_y_ini + tmp_r * sin(tmp_theta) - center_y_rnd) .* cos(tmp_rot) + center_y_rnd;
                
                % computing NND
                D1_ini = pdist2([tmp_x_ini, tmp_y_ini], [CoorFile_1(dumidx, 2), CoorFile_1(dumidx, 3)]);
                NND1_ini = min(D1_ini(:));
                
                % computing randomized NND
                D1_rnd = pdist2([tmp_x_rnd, tmp_y_rnd], [CoorFile_1(dumidx, 2), CoorFile_1(dumidx, 3)]);
                NND1_rnd = min(D1_rnd(:));
                
                fprintf(fid_NND_2, '%f %f %f %f\r\n', ii, numMol, NND1_ini, NND1_rnd);
            end
            fclose(fid_NND_2);
            
            toc;
        end
    end
end
