path0 = uigetdir('', 'Choose the Output Directory');
SampleList = dir(path0);
isSampleDir = [SampleList(:).isdir];
SampleNames = {SampleList(isSampleDir).name}';
SampleNames(ismember(SampleNames, {'.','..','map','Analysis','preview'})) = [];

RoiNum = fopen([path0 '\NumOfROI.txt'],'w');

for i = 1 : numel(SampleNames)
    numROI = 0;
    spools = dir([path0 '\' SampleNames{i} '\*_roi.zip']);
    for j = 1 : numel(spools)
        roi_name = spools(j).name;
        roi_zip = [path0 '\' SampleNames{i} '\' roi_name];
        
        [pth, tempdir, ~] = fileparts(roi_zip);
        tmpROIdir = [pth '\' tempdir];

        if exist(tmpROIdir,'dir') == 7
            rmdir(tmpROIdir, 's');
            mkdir(tmpROIdir);
        else
            mkdir(tmpROIdir);
        end
        unzip(roi_zip, tmpROIdir);

        ROIs = dir([tmpROIdir '\*.roi']);
        numROI = numROI + numel(ROIs);
        rmdir(tmpROIdir, 's');
    end
    fprintf(RoiNum, '%s %d\r\n', SampleNames{i}, numROI);
end
        