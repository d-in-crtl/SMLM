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
imszx = 1200;
imszy = 1200;
CameraPixel = 65.39; % unit = nm, Camera pixel size
fname = 'spool'; % this was set in micromanager during acquisition


%% ========================== Color Information ===========================
AcqSeq = 'DRGB';  % Aquisition Sequence         
NumOfColor = length(AcqSeq);


%% ====================== Rendering Parameters=============================
% Rendering preference:
%   'n' render image with number of fitted emitters within a Target Pixel
%   'b' render image as 1 or 0 if the emitter once or never appeared within a Target Pixel
%   'I' render image with the average intensity of fitted emitters within a Target Pixel
render_mode = 'b'; % n is set to default and recommended                   
TargetPixel = 10; % nm, Image will be reconstructed into images of TargetPixelx pixel size in x- and y- direction


%% ======================= Launch the Fitting =============================
Path0 = uigetdir('', 'Choose the  Output  Directory');

SampleList = dir(Path0); 
isSampleDir = [SampleList(:).isdir];
SampleName = {SampleList(isSampleDir).name}';
SampleName(ismember(SampleName, {'.','..','map','Analysis','preview'})) = [];

% Direct program to read the correct samplefolder sequentially 
for s = 1 : numel(SampleName)    

    SpoolList = dir([Path0 '\' SampleName{s} '\' fname '_*.result']);
    SpoolNames = {SpoolList(:).name}';
    if mod(numel(SpoolNames), NumOfColor) ~= 0
        error('color number and result tables do not matech');
    end
    
    for sp = 1 : numel(SpoolNames) / NumOfColor
        
        if NumOfColor > 1
            out = zeros(ceil(imszy*CameraPixel/TargetPixel), ceil(imszx*CameraPixel/TargetPixel), max(NumOfColor, 3), 'uint16');
        else
            out = zeros(ceil(imszx*CameraPixel/TargetPixel), ceil(imszx*CameraPixel/TargetPixel), 'uint16');
        end
        
        for i = 1 : NumOfColor
            
            fprintf('\n%s\n', ['working on ' SampleName{s} ' ' fname '_' num2str(sp) '_' colordict(AcqSeq(i))]);
            FullFileName = [Path0 '\' SampleName{s} '\' fname '_' num2str(sp) '_' colordict(AcqSeq(i)) '.result'];
            reports = dlmread(FullFileName);
            recon = Coor2Pixel(reports, CameraPixel, TargetPixel, imszx, imszy, render_mode);
            
            out(:, :, colordisp(AcqSeq(i))) = recon;
                
        end
        
        SRname = [Path0 '\' SampleName{s} '\' fname '_' num2str(sp) '_Reconstruction_' num2str(TargetPixel) '_mode_' render_mode '.tif'];
        
        if NumOfColor > 3
            
            t = Tiff(SRname,'w');
            setTag(t, 'ImageLength', ceil(imszx*CameraPixel/TargetPixel));
            setTag(t, 'ImageWidth', ceil(imszx*CameraPixel/TargetPixel));
            setTag(t, 'Photometric', Tiff.Photometric.RGB)
            setTag(t, 'BitsPerSample', 16);
            setTag(t, 'SamplesPerPixel', NumOfColor);
            setTag(t, 'ExtraSamples', 2);
            setTag(t, 'PlanarConfiguration', 1);
            setTag(t, 'Software', 'MATLAB');
            write(t, out);
            close(t);
        else
            imwrite(out, SRname);
        end
        
    end
end

clearvars;
fclose('all');
clear mex;