function [protim3,Dbothfin,segdna,segcell,dnaim3,cellim3,imgcent_coordinate] = getrealimage_hela(imnum)

%{
load(['proc/cell_' num2str(imnum) '/cell' num2str(imnum) '.mat'],'protim3','cellim3','dnaim3');
load(['proc/cell_' num2str(imnum) '/final_info.mat'],'Dbothfin','imgcent_coordinate')
load(['proc/cell_' num2str(imnum) '/man_seg.mat'],'segdna','segcell');
%}
load(['./proc/cell_' num2str(imnum) '/cell' num2str(imnum) '.mat'],'protim3','cellim3','dnaim3');
load(['./proc/cell_' num2str(imnum) '/final_info.mat'],'Dbothfin','imgcent_coordinate')
load(['./proc/cell_' num2str(imnum) '/man_seg.mat'],'segdna','segcell');

a= size(protim3,3);
Dbothfin(:,:,[1:20,a+21:end]) = [];