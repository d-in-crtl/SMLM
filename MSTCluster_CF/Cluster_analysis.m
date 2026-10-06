function [ellipse, S_Ellipse, NumSpots, density] = Cluster_analysis(sub_x, sub_y)

% sub_x and sub_y are coordinates of one cluster. sub_x and sub_y should be in column format
% col_1 and col_2 can be plotted as x and y to draw an ellipse for this cluster

    NumSpots = size(sub_x, 1);
    cluster_center = [mean(sub_x), mean(sub_y)];
    
    % PCA analysis of the cluster
    sub_x = sub_x - cluster_center(1); % adjust the center to 0
    sub_y = sub_y - cluster_center(2);
    cov_M = cov(sub_x, sub_y); % Calculating the covariance matrix
    [V, D] = eig(cov_M); % D is the eigenvalues of cov_M, V is the eigenvectors
    [~, I] = max(max(D)); % first principle component
    
    % Rotating the cluster so that its long arm falls onto the x-axis
    Rot_M = [V(1, I), V(2, I); -V(2, I), V(1, I)]; 
    sub_new = Rot_M * [sub_x'; sub_y'];
    
    % Quantify the rotated cluster (approximate to an ellipse)
    lArm = 2.5*std(sub_new(1, :));
    sArm = 2.5*std(sub_new(2, :));
    S_Ellipse = pi * lArm * sArm;
    density = NumSpots / S_Ellipse;
    
    % ploting the Ellipse and rotate it back
    t = 0:pi/19:2*pi;
    ellipse = Rot_M' * [lArm .* cos(t); sArm .* sin(t)];
    ellipse(1, :) = ellipse(1, :) + cluster_center(1);
    ellipse(2, :) = ellipse(2, :) + cluster_center(2);
    ellipse = ellipse';
    
end

