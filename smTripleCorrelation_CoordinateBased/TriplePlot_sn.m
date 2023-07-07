function [h] = TriplePlot_sn(triple_trans, dim, rho_res, results, coordinates, mask, magnification, TriCalMode, TriDispMode, ColorCode, BlackColor, BlackAlpha, ColorAlpha)

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

% ========================= Setting the dim code ==========================
dim1 = [TriCalMode(1) TriCalMode(3)]; % dist 1_3 indx_1
dim2 = [TriCalMode(1) TriCalMode(2)]; % dist 1_2 indx_2
dim3 = [TriCalMode(2) TriCalMode(3)]; % dist 2_3 indx_3

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

SizeofCata_null = size(coordinates_cata_null, 1);
SizeofCata_0 = size(coordinates_cata_0, 1);
SizeofCata_1 = size(coordinates_cata_1, 1);
SizeofCata_2 = size(coordinates_cata_2, 1);
SizeofCata_3 = size(coordinates_cata_3, 1);

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

% ============================ Plotting ===================================
h = figure('Position', scrsz);
% ---------------- Ploting the heat map into upper pannel ----------------- 
rho_center = ((1:dim)'-0.5) .* rho_res;

TripleCorrelation_12_23 = squeeze(nanmean(triple_trans, 1));
TripleCorrelation_13_23 = squeeze(nanmean(triple_trans, 2));
TripleCorrelation_13_12 = squeeze(nanmean(triple_trans, 3));

subplot(2,12,1:4);
image(rho_center, rho_center, TripleCorrelation_12_23, 'CDataMapping', 'scaled');
axis xy;
colorbar;
set(gca,'FontSize',8,'DataAspectRatio',[1 1 1]);
ylabel(['Distance of ' dim2 ' (nm)'],'FontSize',9);
xlabel(['Distance of ' dim3 ' (nm)'],'FontSize',9);

hold on
plot(results_catanull(:,3), results_catanull(:,2), 'kx', 'MarkerFaceColor', 'k');
plot(results_cata0(:,3), results_cata0(:,2), 'r+', 'MarkerFaceColor', 'r');
plot(results_cata1(:,3), results_cata1(:,2), 'r+', 'MarkerFaceColor', 'r');
plot(results_cata2(:,3), results_cata2(:,2), 'r+', 'MarkerFaceColor', 'r');
plot(results_cata3(:,3), results_cata3(:,2), 'r+', 'MarkerFaceColor', 'r');
hold off
        
subplot(2,12,5:8);
image(rho_center, rho_center, TripleCorrelation_13_23, 'CDataMapping', 'scaled');
axis xy;
colorbar;
set(gca,'FontSize',8,'DataAspectRatio',[1 1 1]);
ylabel(['Distance of ' dim1 ' (nm)'],'FontSize',9);
xlabel(['Distance of ' dim3 ' (nm)'],'FontSize',9);
hold on
plot(results_catanull(:,3), results_catanull(:,1), 'kx', 'MarkerFaceColor', 'k');
plot(results_cata0(:,3), results_cata0(:,1), 'r+', 'MarkerFaceColor', 'r');
plot(results_cata1(:,3), results_cata1(:,1), 'r+', 'MarkerFaceColor', 'r');
plot(results_cata2(:,3), results_cata2(:,1), 'r+', 'MarkerFaceColor', 'r');
plot(results_cata3(:,3), results_cata3(:,1), 'r+', 'MarkerFaceColor', 'r');
hold off
        
subplot(2,12,9:12);
image(rho_center, rho_center, TripleCorrelation_13_12, 'CDataMapping', 'scaled');
axis xy;
colorbar;
set(gca,'FontSize',8,'DataAspectRatio',[1 1 1]);
ylabel(['Distance of ' dim1 ' (nm)'],'FontSize',9);
xlabel(['Distance of ' dim2 ' (nm)'],'FontSize',9);
hold on
plot(results_catanull(:,2), results_catanull(:,1), 'kx', 'MarkerFaceColor', 'k');
plot(results_cata0(:,2), results_cata0(:,1), 'r+', 'MarkerFaceColor', 'r');
plot(results_cata1(:,2), results_cata1(:,1), 'r+', 'MarkerFaceColor', 'r');
plot(results_cata2(:,2), results_cata2(:,1), 'r+', 'MarkerFaceColor', 'r');
plot(results_cata3(:,2), results_cata3(:,1), 'r+', 'MarkerFaceColor', 'r');
hold off

% --------------- Ploting the triangles into lower pannel -----------------
subplot(2,12,13:15);
for k = 1 : SizeofCata_null
    polarplot(phi_cata_null(k, :)', rho_cata_null(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
polarscatter(reshape(phi_cata_null(:, 1:3), [], 1), reshape(rho_cata_null(:, 1:3), [], 1), magnification.*repmat(results_catanull(:,4), 3, 1), 'MarkerFaceColor', BlackColor, 'MarkerFaceAlpha', BlackAlpha, 'MarkerEdgeColor', 'none');

for k = 1 : SizeofCata_0
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

subplot(2,12,16:18);
for k = 1 : SizeofCata_null
    polarplot(phi_cata_null(k, :)', rho_cata_null(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
polarscatter(reshape(phi_cata_null(:, 1:3), [], 1), reshape(rho_cata_null(:, 1:3), [], 1), magnification.*repmat(results_catanull(:,4), 3, 1), 'MarkerFaceColor', BlackColor, 'MarkerFaceAlpha', BlackAlpha, 'MarkerEdgeColor', 'none');

for k = 1 : SizeofCata_1
    polarplot(phi_cata_1(k, :)', rho_cata_1(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_1(:, k), rho_cata_1(:, k), magnification.*results_cata1(:,4), 'MarkerFaceColor', DispCMap(k,:), 'MarkerFaceAlpha', ColorAlpha, 'MarkerEdgeColor', 'none');
end
hold off

subplot(2,12,19:21);
for k = 1 : SizeofCata_null
    polarplot(phi_cata_null(k, :)', rho_cata_null(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
polarscatter(reshape(phi_cata_null(:, 1:3), [], 1), reshape(rho_cata_null(:, 1:3), [], 1), magnification.*repmat(results_catanull(:,4), 3, 1), 'MarkerFaceColor', BlackColor, 'MarkerFaceAlpha', BlackAlpha, 'MarkerEdgeColor', 'none');

for k = 1 : SizeofCata_2
    polarplot(phi_cata_2(k, :)', rho_cata_2(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
for k = 1 : 3
    polarscatter(phi_cata_2(:, k), rho_cata_2(:, k), magnification.*results_cata2(:,4), 'MarkerFaceColor', DispCMap(k,:), 'MarkerFaceAlpha', ColorAlpha, 'MarkerEdgeColor', 'none');
end
hold off

subplot(2,12,22:24);
for k = 1 : SizeofCata_null
    polarplot(phi_cata_null(k, :)', rho_cata_null(k, :)', 'Color', [BlackColor, BlackAlpha]);
    rlim([0 500]);
    rticklabels([]);
    thetaticklabels([]);
    hold on
end
polarscatter(reshape(phi_cata_null(:, 1:3), [], 1), reshape(rho_cata_null(:, 1:3), [], 1), magnification.*repmat(results_catanull(:,4), 3, 1), 'MarkerFaceColor', BlackColor, 'MarkerFaceAlpha', BlackAlpha, 'MarkerEdgeColor', 'none');

for k = 1 : SizeofCata_3
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