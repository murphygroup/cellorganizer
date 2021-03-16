function [randlengths] = getrandlengths_hela(n,mu_len,sigma_len,colli_min_number,cellnum,batchno,subfolder)

%subfolder = 1;

myfile =  ['outputs_' num2str(subfolder) '/images/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/new_n_' num2str(n) '/new_sim_n_' num2str(n) '_mulen_' num2str(mu_len) '_siglen_' num2str(sigma_len) '_colli_' num2str(colli_min_number) '.mat'];

load(myfile,'randlengths');

