function [coors, Cluster_List, Cluster_Figure] = mst_dbscan(points_in, eps, minpts)
    
    % -------------------- MST based DBSCAN clustering --------------------
    % Step 1. Delaunay Triangulation
    DT = delaunayTriangulation(points_in);

    Idx_1 = [DT.ConnectivityList(:, 1); DT.ConnectivityList(:, 1); DT.ConnectivityList(:, 2)];
    Idx_2 = [DT.ConnectivityList(:, 2); DT.ConnectivityList(:, 3); DT.ConnectivityList(:, 3)];
    Idx = sort([Idx_1, Idx_2], 2);
    Idx = unique(Idx, 'rows');

    TriEdges = zeros(size(Idx, 1), 1);
    for ii = 1 : size(Idx, 1)
        TriEdges(ii) = pdist2(DT.Points(Idx(ii, 1), :), DT.Points(Idx(ii, 2), :));
    end

    % Step 2. Searching MinSpanTree from DT edges
    G = graph(Idx(:, 1), Idx(:, 2), TriEdges);
    MST = minspantree(G);

    % Step 3. Cutting the MST at given eps
    mask_1 = MST.Edges.Weight > eps;
    MST = rmedge(MST, find(mask_1));
    bins = conncomp(MST);

    classifications = bins';
    
    % ---------------- Cluster statistics ------------------------
    NumCluster = max(classifications);
    Cluster_List = zeros(NumCluster, 6);
    mask = false(NumCluster, 1);
    
    scrsz = get(groot,'ScreenSize');
    Cluster_Figure = figure('Position', scrsz);
    
    plot(DT.Points(:, 1), DT.Points(:, 2), 'k.');
    hold on
    for cl = 1 : NumCluster
    
        % Extracting coordinates of the cluster
        sub_idx = classifications == cl;
        sub_x = DT.Points(sub_idx, 1);
        sub_y = DT.Points(sub_idx, 2);
        if size(sub_x, 1) >= minpts
            plot(sub_x, sub_y, '.');
        
            % individual cluster analysis
            [ellipse, S_Ellipse, NumSpots, density] = Cluster_analysis(sub_x, sub_y);
            Cluster_List(cl, :) = [cl, NumSpots, S_Ellipse, density, mean(sub_x), mean(sub_y)];
            plot(ellipse(:, 1), ellipse(:, 2), 'r-'); % plotting the ellipse for each cluster
            mask(cl) = true;
        end
    end
    Cluster_List(~mask, :) = [];
    coors = [classifications, DT.Points];
    ax = gca;
    set(ax, 'DataAspectRatioMode', 'Manual', 'DataAspectRatio', [1, 1, 1]);
    
end


