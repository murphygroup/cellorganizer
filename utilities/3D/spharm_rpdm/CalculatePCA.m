function [scales,coeff,score,latent,tsquared,explained,mu,train_score,train_explained,train_coeff] = CalculatePCA(X, latent_dim)
%   
%   Created by Hanxi Xiao & Serena Abraham
%   Copy and Pasted from the Original CellOrganizer
%   train_spharm_rpdm_model.m script


scales = sqrt(sum(X .^ 2, 2));
X = X ./ scales;

[coeff,score,latent,tsquared,explained, mu] = pca(X);

if size(score, 2) < latent_dim
    warning('The given latent dimension %d is larger than the maximum latent dimension in PCA, set latent dimension as the maximum latent dimension %d', latent_dim, size(score, 2));
    latent_dim = size(score, 2);

end

train_score = score(:, 1 : latent_dim);
train_explained = sum(explained(1 : latent_dim));
train_coeff = coeff(:, 1 : latent_dim);

end

