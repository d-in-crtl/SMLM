%% =========================== Initializing ===============================
clearvars;
clc;
fclose('all');
clear mex;

p = mfilename('fullpath');
[filepath, name, ~] = fileparts(p);
cd(filepath);
colordict = containers.Map({'R', 'G', 'B', 'D'}, {'RED', 'GREEN', 'BLUE', 'DarkRed'});
colordisp = containers.Map({'R', 'G', 'B', 'D'}, {1, 2, 3, 4});


%% ======================= Camera Information =============================

CAMERA = 'C:\STORM_A';
DATE = '20190909';

offset = dlmread([CAMERA '\' DATE '\AVG.txt']);
var = dlmread([CAMERA '\' DATE '\VAR.txt']);
gain = dlmread([CAMERA '\' DATE '\GAIN.txt']);

%{
offset = dlmread('C:\Users\rothenberglab\Documents\YY\MATLAB\Camera_B\AVG.txt');
var = dlmread('C:\Users\rothenberglab\Documents\YY\MATLAB\Camera_B\VAR.txt');
gain = dlmread('C:\Users\rothenberglab\Documents\YY\MATLAB\Camera_B\GAIN.txt');
%}

imszx = 1200; % unit = pixel, image width, 1200 by defalt
imszy = 1200; % unit = pixel, image height, 1200 by defalt
CameraPixel = 65.39; % unit = nm, Camera pixel size
fname = 'spool'; % this was set in micromanager during acquisition


%% ========================== Color Information ===========================
AcqSeq = 'R';  % Aquisition Sequence         
RefColor = 'R';  % The color that other colors are mapped to. Must match the RefColor as in channel color mapping step.  

% filtering method for each color channel 
% options: 'wvlet'(default and recommended), 'gauss', 'box'(least preferred)
filtmeth_D = 'wvlet';
filtmeth_R = 'wvlet';
filtmeth_G = 'wvlet';
filtmeth_B = 'wvlet';

% thresholding method for each color channel 
% options: 'std'(default and recommended), a number, times of std of pixel read-out noise
thresh_D = 'std';
thresh_R = 'std';
thresh_G = 'std';
thresh_B = 'std';

NumOfColor = length(AcqSeq);
filtmeth_set = cell(NumOfColor, 1);
thresh_set = cell(NumOfColor, 1);
for i = 1 : NumOfColor
    switch AcqSeq(i)
        case 'D'
            filtmeth_set{i} = filtmeth_D;
            thresh_set{i} = thresh_D;
        case 'R'
            filtmeth_set{i} = filtmeth_R;
            thresh_set{i} = thresh_R;
        case 'G'
            filtmeth_set{i} = filtmeth_G;
            thresh_set{i} = thresh_G;
        case 'B'
            filtmeth_set{i} = filtmeth_B;
            thresh_set{i} = thresh_B;
    end
end


%% ==================== PSF Fitting Parameters=============================
PSFsigma = 145; % nm, expected PSF sigma (initial guess, doesn't have to be exact)
boxsz = 11; % unit = pixels, length of the edge of the square box used for PSF fitting
iterations = 100; % number of iterations used in LBFGSB-MLE fitting method
XFA = 'MFA'; % SFA or MFA for single-emitter fitting or multi-emitter fitting
Algorithm = 'ABFGSB'; % acceptable algorithms: 'NR'-Newton Raphson; 'LM'-Levenberg Marquardt; 'ABFGSB'-Adpated BFGS with Box constration
option = 2; % 1, 2, or 3 for fixed PSF sigma, free PSF sigma, or independent PSF sigma_x and _y
lim = 5E5; % GPU Capacity

% Rendering preference:
%   'n' render image with number of fitted emitters within a Target Pixel
%   'b' render image as 1 or 0 if the emitter once or never appeared within a Target Pixel
%   'I' render image with the average intensity of fitted emitters within a Target Pixel
render_mode = 'n'; % n is set to default and recommended                   
TargetPixel = 10; % nm, Image will be reconstructed into images of TagetPixel pixel size

PSFsigma = PSFsigma / CameraPixel;
zm = CameraPixel / TargetPixel; % Magnification factor in reconstruction.
switch XFA
    case 'SFA'
        Fitter = @GaussSFA2D_SCMOS_wrapper;
    case 'MFA'
        Fitter = @GaussMFA2D_SCMOS_wrapper;
    otherwise
        error('Please specify SFA or MFA');
end


%% ===================== Color Mapping Information ========================
fprintf('The Reference color is set to %s\n', colordict(RefColor));
P_img = cell(NumOfColor, 1);
Q_img = cell(NumOfColor, 1);
P_coor = cell(NumOfColor, 1);
Q_coor = cell(NumOfColor, 1);
for i = 1 : NumOfColor
    if ~strcmp(AcqSeq(i), RefColor)
        [MapName, PathName] = uigetfile('*.map',['Select the' colordict(AcqSeq(i)) 'map file']);
        FullMapName = [PathName, MapName];
        [mappth, mname, ~] = fileparts(FullMapName);
        FullCoorName = [mappth '\' mname '.coor'];
        PQ_img = dlmread(FullMapName);
        P_img{i} = PQ_img(1:size(PQ_img,1)/2);
        Q_img{i} = PQ_img(size(PQ_img,1)/2+1:size(PQ_img,1));
        PQ_coor = dlmread(FullCoorName);
        P_coor{i} = PQ_coor(1:size(PQ_coor,1)/2);
        Q_coor{i} = PQ_coor(size(PQ_coor,1)/2+1:size(PQ_coor,1));
    end
end


%% ======================= Launch the Fitting =============================
% data saving: Path0 --> samplefolder --> spoolfolder --> spool files
Path0 = uigetdir('', 'Choose the  Main  Directory');
[parentpth, ~, ~] = fileparts(Path0);

SampleList = dir(Path0); 
isSampleDir = [SampleList(:).isdir];
SampleName = {SampleList(isSampleDir).name}';
SampleName(ismember(SampleName, {'.','..','map','Analysis','preview'})) = [];

% Create Output folders based on existence of samplefolders
for s = 1 : numel(SampleName)
    if ~exist([parentpth '\Output_2DGauss_' XFA '_' Algorithm '\' SampleName{s}], 'dir')
        mkdir([parentpth '\Output_2DGauss_' XFA '_' Algorithm '\' SampleName{s}]);
    end
end

% Direct program to read the correct samplefolder sequentially 
for s = 1 : numel(SampleName)    
    outputpth = [parentpth '\Output_2DGauss_' XFA '_' Algorithm '\' SampleName{s} '\'];

    SpoolList = dir([Path0 '\' SampleName{s} '\' fname '_*']);
    isSpoolDir = [SpoolList(:).isdir];
    SpoolNames = {SpoolList(isSpoolDir).name}';

    for sp = 1 : NumOfColor : numel(SpoolNames)
        
        if NumOfColor > 1
            out_maxProj = zeros(imszy, imszx, max(NumOfColor, 3), 'uint16');
            out_SR = zeros(ceil(imszx*CameraPixel/TargetPixel), ceil(imszx*CameraPixel/TargetPixel), max(NumOfColor, 3), 'uint16');
        else
            out_maxProj = zeros(imszy, imszx, 'uint16');
            out_SR = zeros(ceil(imszx*CameraPixel/TargetPixel), ceil(imszx*CameraPixel/TargetPixel), 'uint16');
        end
        
        for i = 1 : NumOfColor
            
            fprintf('\n%s\n', ['working on ' SampleName{s} ' ' fname '_' num2str(sp+i-1) '_' colordict(AcqSeq(i))]);
            
            disp('inquiring file information...')
            FileList = dir([Path0 '\' SampleName{s} '\' fname '_' num2str(sp+i-1) '\' fname '_' num2str(sp+i-1) '*.tif']);
            imgnum = zeros(numel(FileList), 1);
            for f = 1 : numel(FileList)
                FullFileName = [Path0 '\' SampleName{s} '\' fname '_' num2str(sp+i-1) '\' FileList(f).name];
                tiffinfo = imfinfo(FullFileName);
                imgnum(f) = numel(tiffinfo);
            end
            disp(['loading images ... (' num2str(sum(imgnum)) ' frames in total)']);
            
            ims = zeros(imszy, imszx, sum(imgnum), 'uint16');
            frmnum = 0;
            for f = 1 : numel(FileList)
                FullFileName = [Path0 '\' SampleName{s} '\' fname '_' num2str(sp+i-1) '\' FileList(f).name];
                warning('off','all');
                t = Tiff(FullFileName, 'r');
                for j = 1 : imgnum(f)
                    t.setDirectory(j);
                    ims(:, :, frmnum+j) = t.read();
                end
                t.close();
                warning('on','all');
                frmnum = frmnum + imgnum(f);
            end

            disp('segmenting and fitting...');
            tic
            out_im = max(ims, [], 3);
            [reports, ims_res, endfrm] = Fitter(ims, 0, offset, var, gain, PSFsigma, iterations, boxsz, lim, filtmeth_set{i}, thresh_set{i}, Algorithm, option);
            while ~isempty(ims_res)
                [reports_res, ims_res, endfrm] = Fitter(ims_res, endfrm, offset, var, gain, PSFsigma, iterations, boxsz, lim, filtmeth_set{i}, thresh_set{i}, Algorithm, option);
                reports = cat(1, reports, reports_res);
                clear reports_res;
            end
            toc
            save([Path0 '\' SampleName{s} '\' fname '_' num2str(sp+i-1) '\' fname '_' num2str(sp+i-1) '.result'], 'reports', '-ASCII');
            imwrite(out_im, [Path0 '\' SampleName{s} '\' fname '_' num2str(sp+i-1) '\' fname '_' num2str(sp+i-1) '.frm1'],'tif');
            
            if ~strcmp(AcqSeq(i), RefColor)
                out_im = movie_correction(out_im, P_img{i}, Q_img{i});
                reports = coor_correction(reports, P_coor{i}, Q_coor{i});
            end
            recon = Coor2Pixel(reports, CameraPixel, TargetPixel, imszx, imszy, render_mode); 
            
            out_maxProj(:, :, colordisp(AcqSeq(i))) = out_im;
            out_SR(:, :, colordisp(AcqSeq(i))) = recon;
            
            save([outputpth fname '_' num2str(ceil(sp/NumOfColor)) '_' colordict(AcqSeq(i)) '.result'],'reports','-ASCII');
                   
        end
        
        if NumOfColor > 3
            t1 = Tiff([outputpth fname '_' num2str(ceil(sp/NumOfColor)) '_frm1.tif'],'w');
            setTag(t1, 'ImageLength', imszy);
            setTag(t1, 'ImageWidth', imszx);
            setTag(t1, 'Photometric', Tiff.Photometric.RGB)
            setTag(t1, 'BitsPerSample', 16);
            setTag(t1, 'SamplesPerPixel', NumOfColor);
            setTag(t1, 'ExtraSamples', 2);
            setTag(t1, 'PlanarConfiguration',1);
            setTag(t1, 'Software','MATLAB');
            write(t1, out_maxProj);
            close(t1);
            
            t2 = Tiff([outputpth fname '_' num2str(ceil(sp/NumOfColor)) '_Reconstruction.tif'],'w');
            setTag(t2, 'ImageLength', ceil(imszx*CameraPixel/TargetPixel));
            setTag(t2, 'ImageWidth', ceil(imszx*CameraPixel/TargetPixel));
            setTag(t2, 'Photometric', Tiff.Photometric.RGB)
            setTag(t2, 'BitsPerSample', 16);
            setTag(t2, 'SamplesPerPixel', NumOfColor);
            setTag(t2, 'ExtraSamples', 2);
            setTag(t2, 'PlanarConfiguration', 1);
            setTag(t2, 'Software', 'MATLAB');
            write(t2, out_SR);
            close(t2);
        else
            imwrite(out_maxProj, [outputpth fname '_' num2str(ceil(sp/NumOfColor)) '_frm1.tif']);
            imwrite(out_SR, [outputpth fname '_' num2str(ceil(sp/NumOfColor)) '_Reconstruction.tif']);
        end
        
    end
end

clearvars;
fclose('all');
clear mex;