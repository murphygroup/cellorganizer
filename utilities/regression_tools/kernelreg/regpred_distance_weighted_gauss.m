function [y_pred, k_x2y, save_all_file] = regpred_distance_weighted_gauss(x, y, cvinds,  savedir, param)
%does nested hold one out cross validation

if ~exist('cvinds', 'var')
    cvinds = 1:size(x,1);
end

if ~exist(savedir, 'dir')
    mkdir(savedir)
end

if ~exist('param', 'var')
    param = [];
end

param = ml_initparam(param, struct('load_dat', true, ...
                                    'build', false));

save_all_file = [savedir filesep 'all_dat.mat'];

if ~exist(save_all_file, 'file') && ~param.build

    kernel_grid = 10.^(-5:0.05:2);

    cellinds = find((~any(isnan(x),2))&(~any(isnan(y),2)));

    ndat = length(cellinds);

    x = x(cellinds,:);
    y = y(cellinds,:);

    dmat_x = squareform(pdist(x));

    [ucvinds, ~, ~] = unique(cvinds);
    cvinds_per_cv = cell(1, length(ucvinds));

    for i = 1:length(ucvinds)
        cvinds_dmat = ucvinds(i) == cvinds;
        dmat_x(cvinds_dmat, cvinds_dmat) = inf;
        
        cvinds_per_cv{i} = find(cvinds_dmat);
    end

    y_pred = nan(size(y, 2));

    for i = 1:ndat
        k_save = [savedir filesep 'k_' num2str(i) '.mat'];
        k_save_tmp = [k_save '.tmp'];
        if ~exist(k_save, 'file') && ~exist(k_save_tmp, 'file')
            system(['touch "' k_save_tmp '"']);

            disp(['pred:' num2str(i) filesep num2str(ndat)])

            traininds = 1:ndat;
            traininds(i) = [];

            k_x2y = find_kernel(dmat_x(traininds, traininds), y(traininds,:), kernel_grid);

            if ~isnan(k_x2y) 
                y_pred = kernel_pred(y, dmat_x(i,:),k_x2y);
            end

            save(k_save, 'y_pred', 'k_x2y')
            
            system(['rm ' k_save_tmp]);
        end
    end
    
    files = dir([savedir filesep 'k_*.mat']);
    
    any_errors = false;
    if length(files) == ndat
        for i = 1:ndat
            k_save = [savedir filesep 'k_' num2str(i) '.mat'];
            try
                load(k_save)
                k_x2y_all(i) = k_x2y;
                y_pred_all(i,:) = y_pred;
            catch
                any_errors = true;
                system(['rm ' k_save])
            end
        end
        
        if ~any_errors
            k_x2y = k_x2y_all;
            y_pred = y_pred_all;

            save(save_all_file, 'y_pred', 'k_x2y')
        end
    end
end

if param.load_dat && exist(save_all_file, 'file')
    load(save_all_file)
else
    y_pred = [];
    k_x2y = [];
end

end

function [k_best, pred_err_sq] = find_kernel(dmat_embed,p2, kernel_grid)

    
    fvals = nan(1, length(kernel_grid));


    for i = 1:length(kernel_grid(:))

        fval = opt_func(p2, dmat_embed, kernel_grid(i));
        fvals(i) = fval;

    end

    [~, ind] = min(fvals);

%     disp('Computing kernel');
    [k_best, fval, exitflag] = fminsearch(@(x) opt_func(p2, dmat_embed, x), kernel_grid(ind), optimset('Display','off', 'MaxIter', 200));

    [~, pred_pos] = opt_func(p2, dmat_embed, k_best);

    pred_err_sq = sum((p2 - pred_pos).^2,2);


    
end

function plot_kernel_error(savename, fvals, kernel_grid)
    
    if ~exist(savename, 'file')
        
        figure('color', 'w'), 
        hold on

        plot(log10(kernel_grid), fvals)

        xlabel('log10(kernel size)')
        ylabel('HOO MSE')

        saveas(gcf, savename, 'png')
        close(gcf)
    end
end

function [diff, pred_pos] = opt_func(pos,dists,h)

    %if not positive definite, the kernel is no good
    if h <= 0
        diff = inf;
        pred_pos = inf;
        return
    end

    pred_pos = kernel_pred(pos, dists, h);
    
    diff = sum(sum((pos-pred_pos).^2,2),1);
    
    diff(isnan(diff)) = inf;

end

function pred_pos = kernel_pred(pos, dists,h)
    w = zeros(size(dists));
    w(:) = mvnpdf(dists(:), 0, h);
%     w(all(w == 0, 2),:) = 1;
    w(all(w==0,2),:) = nan;

    w = w./repmat(sum(w,2), [1, size(w,2)]);
    
    pred_pos = sum(repmat(permute(pos, [3,1,2]), [size(dists,1),1,1]).*repmat(w, [1,1,size(pos,2)]),2);
    pred_pos = permute(pred_pos, [1,3,2]);
    
end

