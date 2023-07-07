clearvars
clc
addpath('Z:\rothenberglab\homes\yiny02\Codes\smTriCorr\CoordinatesBased');
addpath('C:\Users\ERothenberg\Documents\MATLAB\smTriCorr\FFTBased');

Path0 = 'Z:\rothenberglab\homes\yiny02\Codes\smTriCorr\CoordinatesBased\sz_10000\';
corrsz = [1E2, 1E3, 2E3, 5E3, 1E4];

roiszx = 10240;
roiszy = 10240;

benchLog = fopen([Path0 'benchmark.txt'], 'w');
fprintf(benchLog, '%20s %s %s %s %s\r\n', 'density /um^2', 'FFT', 'smGPU', 'smCPP', 'smMatLab');

for sz = 1 : 5
    density = corrsz(sz)/(roiszx*roiszy)*1E6;
    disp(['Working on density ' num2str(density) '/um^2']);
    filepath = [Path0 'sz_' num2str(corrsz(sz)) '\'];
    
    coor_1 = dlmread([filepath 'coor_1.txt']);
    coor_2 = dlmread([filepath 'coor_2.txt']);
    coor_3 = dlmread([filepath 'coor_3.txt']);

    roi = imread([filepath 'roi.tif']);
    img_1 = roi(:, :, 1);
    img_2 = roi(:, :, 2);
    img_3 = roi(:, :, 3);

    tic;
    %[tripleCorr_12, tripleCorr_13, tripleCorr_r_23, tripleCorr_r_32] = ftTripleCorrelation_MatLab(img_1,img_2,img_3, 100);
    FFTtimer = toc;
    %disp(['FFT base method costs ' num2str(FFTtimer) 's']);
    %save([filepath 'FFT_triple_12.txt'],'tripleCorr_12', '-ASCII');
    %save([filepath 'FFT_triple_13.txt'],'tripleCorr_13', '-ASCII');
    %save([filepath 'FFT_triple_23.txt'],'tripleCorr_r_23', '-ASCII');
    %save([filepath 'FFT_triple_32.txt'],'tripleCorr_r_32', '-ASCII');
    
    tic;
    [Correlation_GPU] = smTripleCorrelation_GPU(coor_1, coor_2, coor_3, roiszx, roiszy);
    GPUtimer = toc;
    disp(['Coordinates base method (GPU) costs ' num2str(GPUtimer) 's']);
    Correlation_GPU = reshape(Correlation_GPU, 90, 100, 100);
    Correlation_GPU = fftshift(Correlation_GPU, 1);
    save([filepath 'Correlation_GPU.mat'], 'Correlation_GPU');
    clear mex;
    
    tic;
    [Correlation_CPP] = smTripleCorrelation_CPP(coor_1, coor_2, coor_3, roiszx, roiszy);
    CPPtimer = toc;
    disp(['Coordinates base method (CPU_C++) costs ' num2str(CPPtimer) 's']);
    Correlation_CPP = reshape(Correlation_CPP, 90, 100, 100);
    Correlation_CPP = fftshift(Correlation_CPP, 1);
    save([filepath 'Correlation_CPP.mat'], 'Correlation_CPP');
    clear mex;
    
    tic;
    [Correlation_MatLab, rho_center, rad_center] = smTripleCorrelation_MatLab(coor_1, coor_2, coor_3, roiszx, roiszy, 100, 5, 2*pi/90);
    MatLabtimer = toc;
    disp(['Coordinates base method (CPU_Matlab) costs ' num2str(MatLabtimer) 's']);
    save([filepath 'Correlation_Matlab.mat'], 'Correlation_MatLab');
    
    fprintf(benchLog, '%f %f %f %f %f\r\n', corrsz(sz), FFTtimer, GPUtimer, CPPtimer, MatLabtimer);
end
fclose(benchLog);