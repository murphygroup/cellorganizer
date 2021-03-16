function [finTable] = recover_parameters_from_real_image_ownlib(subfolder,batchno,cellnum,w)

imCategory = 'general';
featType = 'all';


load(['outputs_' num2str(subfolder) '/featvals/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/' imCategory '/' featType '/forpart2_all_feats_G_psf_batch_' num2str(batchno) '.mat'],'fin_mat');
feat_vector_MI = fin_mat(:,[5:end]);

%subfolder= 1;
%w = [1;1;1;1;1;1;1];
distmet = 'nvar';

[protim3,Dbothfin,segdna,segcell,dnaim3,cellim3,imgcent_coordinate] = getrealimage_hela(cellnum);
feat_vector_ori = real_im_feat_extract(protim3,imgcent_coordinate,w,cellnum);

load allrealfeats all_feat_vector2 idxes

ffidx = [];
for I = 1:length(w)
	if w(I) == 1
		ffidx = [ffidx,idxes{I}(1):idxes{I}(2)];
	end
end

feat_vector = all_feat_vector2(:,ffidx);
feat_vector_MI = feat_vector_MI(:,ffidx);  %%
idxremove =  find(sum(abs(feat_vector),1)==0);

% compute distances and pick the global minimum

% corrections

feat_vector(:,idxremove) = [];
feat_vector_MI(:,idxremove) = [];
feat_vector_ori(idxremove) = [];


varmat = var(feat_vector);
for T = 1:size(fin_mat,1)
	mahdist(T) = nvardist_extract(feat_vector_MI(T,:),feat_vector_ori,varmat);
end
		
[I,J] = min(mahdist);
		
finTable = [fin_mat(J,1:4), I];

end % End of function
