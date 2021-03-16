function [y_pred, k_x2y, save_all_file] = regpred_linreg(x, y, cvinds, savedir, param)
%does nested hold one out cross validation

    if ~exist('cvinds', 'var') | isempty(cvinds)
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

    keepinds = find((~any(isnan(x),2))&(~any(isnan(y),2)));

    ndat = length(keepinds);

    x = x(keepinds,:);
    y = y(keepinds,:);
    
    k_x2y = cell(ndat, 1);
    
    y_pred = nan(ndat, size(y,2));
    
    save_all_file = [savedir filesep 'all_dat'];

    if ~exist([save_all_file '.mat'], 'file') && ~param.build
        [can_start, ~, final_exists, temp_name] = chunk_start( save_all_file);

        if can_start && ~final_exists

            for i = 1:ndat
%                 disp(['pred:' num2str(i) filesep num2str(ndat)])

                traininds = 1:ndat;
                traininds(cvinds == cvinds(i)) = [];

                try
                    k_x2y{i}= find_kernel(x(traininds,:), y(traininds,:));

                    y_pred(i,:) = kernel_pred(x(i,:),k_x2y{i});
                catch
                    k_x2y{i} = nan;
                    y_pred(i,:) = nan;
                end

            end
            save(save_all_file, 'y_pred', 'k_x2y')

            chunk_finish(temp_name)
        end
    end

    save_all_file = [save_all_file '.mat'];

    if exist(save_all_file, 'file') && param.load_dat
        load(save_all_file)
    else
        y_pred = [];
        k_x2y = [];
    end
end

function [beta, pred_err_sq] = find_kernel(x,y)
    x = [ones(size(x,1), 1) x];
    beta = (x'*x)\x'*y;
    
    pred_err_sq = sum((x*beta - y).^2,2);

%     if plot_err
%         plot_kernel_error(savefigurepath, fvals, kernel_grid)
%     end
    
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


function pred_pos = kernel_pred(pos, beta)
    pred_pos = [ones(size(pos,1),1) pos] * beta;
end


