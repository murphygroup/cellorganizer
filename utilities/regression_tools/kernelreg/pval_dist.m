function [ x_pvals, error_norm ] = pval_dist(errors, error_distr )
%CELL_NUC_PRED_PVAL2 Summary of this function goes here
%   Detailed explanation goes here

errors = errors.^2;
error_distr = error_distr.^2;

further_than = repmat(errors, [1, size(error_distr,2)]) >= error_distr;
further_than(logical(eye(size(further_than)))) = 0;

x_pvals = (sum(further_than,2) + 1)./(size(further_than,2)-1 + 1);


error_norm = errors./mean(error_distr,2);


end

