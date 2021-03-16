function [r_lowd_data, r_pca] = tcell_pca( p_data, p_min_var )
% 2011-12-02 tebuck: copied from teb_pca.m.
% 2015-11-02 xruan: change the name of the function to tcell_pca in case of conflection 
% with the built-in pca function. 
min_var = 0.95;
if ( nargin > 1 )
	min_var = p_min_var;
end

p_data = zscore(p_data);
[data_coeff, data_scores, data_latent] = princomp(p_data);

num_pcs = 1;
sum_latent = sum(data_latent);
while(sum(data_latent(1:num_pcs))/sum_latent < min_var)
	num_pcs = num_pcs + 1;
end

%r_lowd_data = data_scores * data_coeff(:, 1:num_pcs);
r_lowd_data = data_scores(:, 1:num_pcs);
%r_lowd_data2 = p_data * data_coeff(:, 1:num_pcs);

%mean(mean((r_lowd_data).^2))
%mean(mean((r_lowd_data2).^2))
%mean(mean((r_lowd_data2 - r_lowd_data).^2))

%clf
%hold on;
%plot( r_lowd_data(:, 1), r_lowd_data(:, 2), 'rx' );
%plot( r_lowd_data2(:, 1), r_lowd_data2(:, 2), 'bo' );
%legend( {'Standard projection', 'My projection'} );
%hold off;

r_pca.coeff = data_coeff;
r_pca.latent = data_latent;
r_pca.num_pcs = num_pcs;
r_pca.pcs_prop_var = sum(data_latent(1:num_pcs))/sum_latent;


	
