function [h_allbut0, h0, h1, h2, h3] = TriplePlot_All(results, coordinates_raw, coordinates_norm, mask, magnification, TriDispMode, ColorCode, BlackColor, BlackAlpha, ColorAlpha, PlotMode)

scrsz = get(groot,'ScreenSize');
% ========================= Setting the color code ========================
DispCMap = zeros(3, 3);
for cc = 1 : 3
    cl = TriDispMode(cc);
    switch cl
        case 'R'
            DispCMap(cc, :) = ColorCode(1, :);
        case 'G'
            DispCMap(cc, :) = ColorCode(2, :);
        case 'B'
            DispCMap(cc, :) = ColorCode(3, :);
        case 'D'
            DispCMap(cc, :) = ColorCode(4, :);
        otherwise
            error(['the color ' cl ' is not coded here']);
    end
end

% ========================== Plotting mode ================================
switch PlotMode
    case 'norm'
        coordinates = coordinates_norm;
    case 'raw'
        coordinates = coordinates_raw;
    otherwise
        coordinates = coordinates_norm;
end

% =========================== Catagorize ==================================
cata_null = mask == -1;
cata_0 = mask == 0;
cata_1 = mask == 1;
cata_2 = mask == 2;
cata_3 = mask == 3;

results_catanull = results(cata_null, :);
results_cata0 = results(cata_0, :);
results_cata1 = results(cata_1, :);
results_cata2 = results(cata_2, :);
results_cata3 = results(cata_3, :);

coordinates_cata_null = coordinates(cata_null, :);
coordinates_cata_0 = coordinates(cata_0, :);
coordinates_cata_1 = coordinates(cata_1, :);
coordinates_cata_2 = coordinates(cata_2, :);
coordinates_cata_3 = coordinates(cata_3, :);

SizeofCata_null = sum(cata_null, 1);
SizeofCata_0 = sum(cata_0, 1);
SizeofCata_1 = sum(cata_1, 1);
SizeofCata_2 = sum(cata_2, 1);
SizeofCata_3 = sum(cata_3, 1);

% ========================== Car2pol transform ============================
x_cata_null = [coordinates_cata_null(:,1), coordinates_cata_null(:,3), coordinates_cata_null(:,5), coordinates_cata_null(:,1)];
y_cata_null = [coordinates_cata_null(:,2), coordinates_cata_null(:,4), coordinates_cata_null(:,6), coordinates_cata_null(:,2)];

x_cata_0 = [coordinates_cata_0(:,1), coordinates_cata_0(:,3), coordinates_cata_0(:,5), coordinates_cata_0(:,1)];
y_cata_0 = [coordinates_cata_0(:,2), coordinates_cata_0(:,4), coordinates_cata_0(:,6), coordinates_cata_0(:,2)];

x_cata_1 = [coordinates_cata_1(:,1), coordinates_cata_1(:,3), coordinates_cata_1(:,5), coordinates_cata_1(:,1)];
y_cata_1 = [coordinates_cata_1(:,2), coordinates_cata_1(:,4), coordinates_cata_1(:,6), coordinates_cata_1(:,2)];

x_cata_2 = [coordinates_cata_2(:,1), coordinates_cata_2(:,3), coordinates_cata_2(:,5), coordinates_cata_2(:,1)];
y_cata_2 = [coordinates_cata_2(:,2), coordinates_cata_2(:,4), coordinates_cata_2(:,6), coordinates_cata_2(:,2)];

x_cata_3 = [coordinates_cata_3(:,1), coordinates_cata_3(:,3), coordinates_cata_3(:,5), coordinates_cata_3(:,1)];
y_cata_3 = [coordinates_cata_3(:,2), coordinates_cata_3(:,4), coordinates_cata_3(:,6), coordinates_cata_3(:,2)];

phi_cata_null = zeros(SizeofCata_null, 4);
rho_cata_null = zeros(SizeofCata_null, 4);

phi_cata_0 = zeros(SizeofCata_0, 4);
rho_cata_0 = zeros(SizeofCata_0, 4);

phi_cata_1 = zeros(SizeofCata_1, 4);
rho_cata_1 = zeros(SizeofCata_1, 4);

phi_cata_2 = zeros(SizeofCata_2, 4);
rho_cata_2 = zeros(SizeofCata_2, 4);

phi_cata_3 = zeros(SizeofCata_3, 4);
rho_cata_3 = zeros(SizeofCata_3, 4);

for k = 1 : 4
    [phi_cata_null(:, k), rho_cata_null(:, k)] = cart2pol(x_cata_null(:, k), y_cata_null(:, k));
    [phi_cata_0(:, k), rho_cata_0(:, k)] = cart2pol(x_cata_0(:, k), y_cata_0(:, k));
    [phi_cata_1(:, k), rho_cata_1(:, k)] = cart2pol(x_cata_1(:, k), y_cata_1(:, k));
    [phi_cata_2(:, k), rho_cata_2(:, k)] = cart2pol(x_cata_2(:, k), y_cata_2(:, k));
    [phi_cata_3(:, k), rho_cata_3(:, k)] = cart2pol(x_cata_3(:, k), y_cata_3(:, k));
end

results_catanull(rho_cata_null(:, 3) > 500, :) = [];
phi_cata_null(rho_cata_null(:, 3) > 500, :) = [];
rho_cata_null(rho_cata_null(:, 3) > 500, :) = [];

results_cata0(rho_cata_0(:, 3) > 500, :) = [];
phi_cata_0(rho_cata_0(:, 3) > 500, :) = [];
rho_cata_0(rho_cata_0(:, 3) > 500, :) = [];

results_cata1(rho_cata_1(:, 3) > 500, :) = [];
phi_cata_1(rho_cata_1(:, 3) > 500, :) = [];
rho_cata_1(rho_cata_1(:, 3) > 500, :) = [];

results_cata2(rho_cata_2(:, 3) > 500, :) = [];
phi_cata_2(rho_cata_2(:, 3) > 500, :) = [];
rho_cata_2(rho_cata_2(:, 3) > 500, :) = [];

results_cata3(rho_cata_3(:, 3) > 500, :) = [];
phi_cata_3(rho_cata_3(:, 3) > 500, :) = [];
rho_cata_3(rho_cata_3(:, 3) > 500, :) = [];

% =================== Plotting all but 0 into h_allbut0 ===================
h_allbut0 = figure('Position', scrsz);

for k = 1 : size(phi_cata_null, 1)
    polarplot(phi_cata_null(k, :)', rho_cata_null(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_null(:, k), rho_cata_null(:, k), magnification.*results_catanull(:,4), 'MarkerFaceColor', DispCMap(k,:), 'MarkerFaceAlpha', ColorAlpha, 'MarkerEdgeColor', 'none');
    %polarscatter(phi_cata_null(:, k), rho_cata_null(:, k), magnification.*log10(results_catanull(:,4)), 'MarkerFaceColor', DispCMap(k,:), 'MarkerFaceAlpha', ColorAlpha, 'MarkerEdgeColor', 'none');
end

for k = 1 : size(phi_cata_1, 1)
    polarplot(phi_cata_1(k, :)', rho_cata_1(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_1(:, k), rho_cata_1(:, k), magnification.*results_cata1(:,4), 'MarkerFaceColor', DispCMap(k,:), 'MarkerFaceAlpha', ColorAlpha, 'MarkerEdgeColor', 'none');
end

for k = 1 : size(phi_cata_2, 1)
    polarplot(phi_cata_2(k, :)', rho_cata_2(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_2(:, k), rho_cata_2(:, k), magnification.*results_cata2(:,4), 'MarkerFaceColor', DispCMap(k,:), 'MarkerFaceAlpha', ColorAlpha, 'MarkerEdgeColor', 'none');
end

for k = 1 : size(phi_cata_3, 1)
    polarplot(phi_cata_3(k, :)', rho_cata_3(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_3(:, k), rho_cata_3(:, k), magnification.*results_cata3(:,4), 'MarkerFaceColor', DispCMap(k,:), 'MarkerFaceAlpha', ColorAlpha, 'MarkerEdgeColor', 'none');
end
hold off

% ========================= Plotting cata_0 into h0 =======================
h0 = figure('Position', scrsz);

for k = 1 : size(phi_cata_0, 1)
    polarplot(phi_cata_0(k, :)', rho_cata_0(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_0(:, k), rho_cata_0(:, k), magnification.*results_cata0(:,4), 'MarkerFaceColor', DispCMap(k,:), 'MarkerFaceAlpha', ColorAlpha, 'MarkerEdgeColor', 'none');
end
hold off

% ========================= Plotting cata_1 into h1 =======================
h1 = figure('Position', scrsz);

for k = 1 : size(phi_cata_null, 1)
    polarplot(phi_cata_null(k, :)', rho_cata_null(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
polarscatter(reshape(phi_cata_null(:, 1:3), [], 1), reshape(rho_cata_null(:, 1:3), [], 1), magnification.*repmat(results_catanull(:,4), 3, 1), 'MarkerFaceColor', BlackColor, 'MarkerFaceAlpha', BlackAlpha, 'MarkerEdgeColor', 'none');

for k = 1 : size(phi_cata_1, 1)
    polarplot(phi_cata_1(k, :)', rho_cata_1(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_1(:, k), rho_cata_1(:, k), magnification.*results_cata1(:,4), 'MarkerFaceColor', DispCMap(k,:), 'MarkerFaceAlpha', ColorAlpha, 'MarkerEdgeColor', 'none');
end

for k = 1 : size(phi_cata_2, 1)
    polarplot(phi_cata_2(k, :)', rho_cata_2(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_2(:, k), rho_cata_2(:, k), magnification.*results_cata2(:,4), 'MarkerFaceColor', BlackColor, 'MarkerFaceAlpha', BlackAlpha, 'MarkerEdgeColor', 'none');
end

for k = 1 : size(phi_cata_3, 1)
    polarplot(phi_cata_3(k, :)', rho_cata_3(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_3(:, k), rho_cata_3(:, k), magnification.*results_cata3(:,4), 'MarkerFaceColor', BlackColor, 'MarkerFaceAlpha', BlackAlpha, 'MarkerEdgeColor', 'none');
end
hold off

% ========================= Plotting cata_2 into h2 =======================
h2 = figure('Position', scrsz);

for k = 1 : size(phi_cata_null, 1)
    polarplot(phi_cata_null(k, :)', rho_cata_null(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
polarscatter(reshape(phi_cata_null(:, 1:3), [], 1), reshape(rho_cata_null(:, 1:3), [], 1), magnification.*repmat(results_catanull(:,4), 3, 1), 'MarkerFaceColor', BlackColor, 'MarkerFaceAlpha', BlackAlpha, 'MarkerEdgeColor', 'none');

for k = 1 : size(phi_cata_1, 1)
    polarplot(phi_cata_1(k, :)', rho_cata_1(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_1(:, k), rho_cata_1(:, k), magnification.*results_cata1(:,4), 'MarkerFaceColor', BlackColor, 'MarkerFaceAlpha', BlackAlpha, 'MarkerEdgeColor', 'none');
end

for k = 1 : size(phi_cata_2, 1)
    polarplot(phi_cata_2(k, :)', rho_cata_2(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_2(:, k), rho_cata_2(:, k), magnification.*results_cata2(:,4), 'MarkerFaceColor', DispCMap(k,:), 'MarkerFaceAlpha', ColorAlpha, 'MarkerEdgeColor', 'none');
end

for k = 1 : size(phi_cata_3, 1)
    polarplot(phi_cata_3(k, :)', rho_cata_3(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_3(:, k), rho_cata_3(:, k), magnification.*results_cata3(:,4), 'MarkerFaceColor', BlackColor, 'MarkerFaceAlpha', BlackAlpha, 'MarkerEdgeColor', 'none');
end
hold off

% ========================= Plotting cata_3 into h3 =======================
h3 = figure('Position', scrsz);

for k = 1 : size(phi_cata_null, 1)
    polarplot(phi_cata_null(k, :)', rho_cata_null(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
polarscatter(reshape(phi_cata_null(:, 1:3), [], 1), reshape(rho_cata_null(:, 1:3), [], 1), magnification.*repmat(results_catanull(:,4), 3, 1), 'MarkerFaceColor', BlackColor, 'MarkerFaceAlpha', BlackAlpha, 'MarkerEdgeColor', 'none');

for k = 1 : size(phi_cata_1, 1)
    polarplot(phi_cata_1(k, :)', rho_cata_1(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_1(:, k), rho_cata_1(:, k), magnification.*results_cata1(:,4), 'MarkerFaceColor', BlackColor, 'MarkerFaceAlpha', BlackAlpha, 'MarkerEdgeColor', 'none');
end

for k = 1 : size(phi_cata_2, 1)
    polarplot(phi_cata_2(k, :)', rho_cata_2(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_2(:, k), rho_cata_2(:, k), magnification.*results_cata2(:,4), 'MarkerFaceColor', BlackColor, 'MarkerFaceAlpha', BlackAlpha, 'MarkerEdgeColor', 'none');
end

for k = 1 : size(phi_cata_3, 1)
    polarplot(phi_cata_3(k, :)', rho_cata_3(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_3(:, k), rho_cata_3(:, k), magnification.*results_cata3(:,4), 'MarkerFaceColor', DispCMap(k,:), 'MarkerFaceAlpha', ColorAlpha, 'MarkerEdgeColor', 'none');
end
hold off