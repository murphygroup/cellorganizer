function [y_pred_out] = kernel_pred_iter(x, y, cvinds, niter, saveparent, param)
%does nested hold one out cross validation

if ~exist('param', 'var')
    param = [];
end
    
param = ml_initparam(param, struct('method', 'gauss', ...
                                    'load_dat', false, ...
                                    'normfunc', @(x) x, ...
                                    'split_on_sig', false, ...
                                    'build', false, ...
                                    'subsample', niter));

results_file = [saveparent filesep 'results.mat'];

if exist(results_file, 'file')
    try
        load(results_file);
        return
    catch
        delete(results_file)
    end
end


    
param.load_dat = false;

ndat = size(x,1);

switch param.method
    case 'gauss'
        [~, ~, save_pred_file] = regpred_distance_weighted_gauss(x, y, cvinds, saveparent, param);
    case 'linreg'
        [~, ~, save_pred_file] = regpred_linreg(x, y, cvinds, saveparent, param);
end

save_perm_file = cell(niter,1);

rng('shuffle')
% iters = randperm(niter);
iters = 1:niter;

%pick a random iteration to work on
for i = 1:length(iters)
    disp(['Iter ' num2str(i) filesep num2str(niter)]);
    save_perm_dir = [saveparent filesep 'perm' num2str(iters(i))];

    %set the random to the iteration index
    rng(iters(i));
    %shuffle both the cell positions and nuc positions cause why not
    x_perm = x(randperm(ndat),:);

    switch param.method
        case 'gauss'
            [~, ~, save_perm_file{iters(i)}] = regpred_distance_weighted_gauss(x_perm, y, cvinds, save_perm_dir, param);
        case 'linreg'
            [~, ~, save_perm_file{iters(i)}] = regpred_linreg(x_perm, y, cvinds, save_perm_dir, param);
    end 
end


% y_m = mean(y(:));
% y_std = std(y(:)-y_m(:));
% y_norm = (y - y_m)./y_std;

if exist(save_pred_file, 'file')
    pred_dat = load(save_pred_file);
else
    y_pred_out = [];
    return
end
y_pred = pred_dat.y_pred;
k_x2y = pred_dat.k_x2y;

y_pred = param.normfunc(y_pred);


% y_pred_norm = (y_pred - y_m)./y_std;
y_pred_err = sum((y - y_pred).^2,2);



y_perm_err = nan(ndat, niter);
for i = 1:niter
    try
        y_pred_perm{i} = load(save_perm_file{i});
        y_iter = y_pred_perm{i}.y_pred;
        
        y_iter = param.normfunc(y_iter);
%         y_iter_norm = (y_iter - y_m)./y_std;

        y_perm_err(:,i) = sum((y - y_iter).^2,2);
        
    catch
        y_perm_err(:,i) = nan;
        
    end
end

keepiter = find(~all(isnan(y_perm_err),1));

if param.subsample < length(keepiter)
    keepiter = keepiter(1:param.subsample);
end
    

isdone = length(keepiter)/niter;





keepdat = ~any(isnan(y_perm_err(:,keepiter)),2);

worse_than_rand = repmat(y_pred_err(keepdat), [1,length(keepiter)]) >= y_perm_err(keepdat,keepiter);


y_perm_err_temp = y_perm_err(keepdat,keepiter);
mean_rand_err = mean(y_perm_err_temp(:));

y_pval = nan(ndat,1);
y_pval(keepdat,:) = (sum(worse_than_rand,2)+eps)/(size(worse_than_rand,2)+eps);

y_pval_all =  (sum(worse_than_rand(:))+eps)/(numel(worse_than_rand(:))+eps);

y_pred_out.x = x;
y_pred_out.y = y;

y_pred_out.y_pred       = y_pred;
y_pred_out.y_pred_err   = y_pred_err;
y_pred_out.y_pval       = y_pval;
y_pred_out.y_pval_all   = y_pval_all;

y_pred_out.y_pval_prod = prod(y_pval);

y_pred_out.k_x2y = k_x2y;

y_pred_out.param        = param;
y_pred_out.isdone       = isdone;

y_pred_out.mean_rand_err = mean_rand_err;

if ~exist('y_pred_perm', 'var')
    y_pred_perm = nan;
end

y_pred_out.y_pred_perm = y_pred_perm;

    



%split the data into statistically significant and not, and attempt to
%rebuild the models
if param.split_on_sig
    param.split_on_sig = false;
    
    stat_sig = y_pred_out.y_pval < 0.05;
    
    if any(stat_sig)
        [y_pred_out_sig] = kernel_pred_iter(x(stat_sig,:), y(stat_sig,:), cvinds(stat_sig), niter, [saveparent '_sig'], param);
    else
        y_pred_out_sig = [];
    end
    
    if any(~stat_sig)
        [y_pred_out_not_sig] = kernel_pred_iter(x(~stat_sig,:), y(~stat_sig,:), cvinds(~stat_sig), niter, [saveparent '_not_sig'], param);
    else
        y_pred_out_not_sig = [];
    end
    
    
    param.split_on_sig = true;
    
    y_pred_out.stat_sig = stat_sig;
    y_pred_out.split_sig = y_pred_out_sig;
    y_pred_out.split_not_sig = y_pred_out_not_sig;
    
end

if isdone == 1
    save(results_file, 'y_pred_out')
end
