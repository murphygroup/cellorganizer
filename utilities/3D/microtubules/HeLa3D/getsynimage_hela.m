function [G_psf,imgcent_coordinate,imXYZ,G,mtXYZ] = getsynimage_hela(n,mu_len,sigma_len,colli_min_number,cellnum,batchno,subfolder,psf)


%myfile =  ['outputs_' num2str(subfolder) '/images/cell_' num2str(cellnum)
%'/batch_' num2str(batchno) '/new_n_' num2str(n) '/new_sim_n_' num2str(n)
%'_mulen_' num2str(mu_len) '_siglen_' num2str(sigma_len) '_colli_'
%num2str(colli_min_number) '.mat'];% DPS 2/15/12 changed path
%myfile =  ['inter_results/HeLa3D/outputs_' num2str(subfolder)
%'/images/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/new_n_'
%num2str(n) '/new_sim_n_' num2str(n) '_mulen_' num2str(mu_len) '_siglen_'
%num2str(sigma_len) '_colli_' num2str(colli_min_number) '.mat']; %DPS
%2/15/12 added '-tcheck' to path
%myfile =  ['inter_results/HeLa3D/outputs_' num2str(subfolder)
%'/images/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/new_n_'
%num2str(n) '/new_sim_n_' num2str(n) '_mulen_' num2str(mu_len) '_siglen_'
%num2str(sigma_len) '_colli_' num2str(colli_min_number) '-tcheck.mat'];%DPS
%DPS 2/15/12 '-tcheck' files do not have full G parameters 


myfile =  ['inter_results/HeLa3D/outputs_' num2str(subfolder) '/images/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/new_n_' num2str(n) '/new_sim_n_' num2str(n) '_mulen_' num2str(mu_len) '_siglen_' num2str(sigma_len) '_colli_' num2str(colli_min_number) '.mat'];

load(myfile,'imgcent_coordinate','randlengths','imXYZ','G','mtXYZ');
G(:,:,[(1:20),(size(G,3)-19:end)]) = [];
%%DPS 2/15/12 added this if statement to allow for non-psf tag
if psf == 1
  G_psf = psf_blur_hela_mean(G);
  G_psf = setsingMTinten_hela(G_psf);
else 
    G_psf = G;
end


imgcent_coordinate(3) = imgcent_coordinate(3) - 20;
