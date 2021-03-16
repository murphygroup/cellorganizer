function [G_psf,imgcent_coordinate,imXYZ,G,mtXYZ] = getsynimage_hela_dps(n,mu_len,colli_min_number,psffile)


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


myfile =  ['./outputs/images/new_n_' num2str(n) '/new_sim_n_' num2str(n) '_mulen_' num2str(mu_len) '_colli_' num2str(colli_min_number) '.mat'];

load(myfile,'imgcent_coordinate','randlengths','imXYZ','G','mtXYZ');
G(:,:,[(1:20),(size(G,3)-19:end)]) = [];
%%DPS 2/15/12 added this if statement to allow for non-psf tag
if size(psffile)>0
  G_psf = psf_blur_hela_mean_dps(G,psffile);
  G_psf = setsingMTinten_hela(G_psf);
else 
    G_psf = G;
end


imgcent_coordinate(3) = imgcent_coordinate(3) - 20;
