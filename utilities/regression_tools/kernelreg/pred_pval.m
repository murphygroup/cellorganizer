function [ x_pvals, norm_error ] = pred_pval( x, x_pred )
%PRED_PVAL Summary of this function goes here
%   Detailed explanation goes here

error_distr = squareform(pdist(x));

errors = sqrt(sum((x - x_pred).^2,2));

[x_pvals, norm_error] = pval_dist(errors, error_distr);

end

