clearvars
clc
fclose('all');

p = mfilename('fullpath');
[filepath, name, ~] = fileparts(p);
addpath(filepath);

load([filepath '\test.mat']);

eps = 50; % Unit = nm
minpts = 10;

disp('Analyzing...');
[classifications, Cluster_List, Cluster_Figure] = mst_dbscan(points, eps, minpts);
disp('Done');

