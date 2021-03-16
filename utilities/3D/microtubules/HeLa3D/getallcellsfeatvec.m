clear all
close all

addpath slicfiles/

%w = [1;1;1;1;1;1;1];
w = [1;1;1;1;1;1;1;1;1;1;1;1;1];

imCategory = 'general';
I = 1;
for cellnum = [1 2 4 6 10 12:42 46:49 51:52]
	[protim3,Dbothfin,segdna,segcell,dnaim3,cellim3,imgcent_coordinate] = getrealimage_hela(cellnum);
 	[all_feat_vector(cellnum,:),idxes] = getfeatvector(protim3,imgcent_coordinate,imCategory,w,cellnum);
	all_feat_vector2(I,:) = all_feat_vector(cellnum,:);
	I = I + 1;
        cellnum
end


save allrealfeats all_feat_vector2 idxes
