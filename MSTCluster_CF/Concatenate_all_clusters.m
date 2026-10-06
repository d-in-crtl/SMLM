clear all

cpath=pwd;
[fname,fpath]=uigetfile('*.Cluster', 'MultiSelect', 'on');
Area_ellipse = [];

for i = 1:size(fname,2)
    data = dlmread(fname{i},' ',1,0);
    Area_ellipse = [Area_ellipse; data(:,3)];
end