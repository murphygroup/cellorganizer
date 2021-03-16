function [ output_args ] = test_reconstruction_methods( distances_complete, savedir )
%Distances_complete may be a sparse incomplete distance matrix

rng(1)

if ~exist('savedir', 'var') | isempty(savedir)
    savedir = pwd;
end

if ~exist(savedir, 'dir')
    mkdir(savedir)
end

tempdir = [savedir filesep 'temp'];
if ~exist(tempdir, 'dir')
    mkdir(tempdir)
end

methods = {struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', inf), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 15), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 14), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 13), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 12), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 11), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 10), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 9), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 8), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 7), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 6), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 5), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 4), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 3), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 2), ...
    struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', 1)};
%, ...;
%             struct('method', 'shortest path',   'weight_factor', 0, 'desired_shape_space_dimensionality', inf), ...
%             struct('method', 'shortest path',   'weight_factor', 1, 'desired_shape_space_dimensionality', inf), ...
%             struct('method', 'shortest path',   'weight_factor', 2, 'desired_shape_space_dimensionality', inf)};%, ...
%
%             struct('method', 'regress',         'weight_factor', 1, 'desired_shape_space_dimensionality', inf), ...
%             struct('method', 'regress',         'weight_factor', 2, 'desired_shape_space_dimensionality', inf), ...
%             struct('method', 'regress',         'weight_factor', 0, 'desired_shape_space_dimensionality', inf), ...
%             struct('method',  'shortest path',   'weight_factor', 0, 'desired_shape_space_dimensionality', 7), ...
%             struct('method', 'shortest path',   'weight_factor', 1, 'desired_shape_space_dimensionality', 7), ...
%             struct('method', 'shortest path',   'weight_factor', 2, 'desired_shape_space_dimensionality', 7), ...
%             struct('method', 'regress',         'weight_factor', 0, 'desired_shape_space_dimensionality', 7), ...
%             struct('method', 'regress',         'weight_factor', 1, 'desired_shape_space_dimensionality', 7), ...
%             struct('method', 'regress',         'weight_factor', 2, 'desired_shape_space_dimensionality', 7), ...
%             };

% methods = {struct('method',  'shortest path',   'weight_factor', 0, 'desired_shape_space_dimensionality', inf), ...
%             struct('method', 'shortest path',   'weight_factor', 1, 'desired_shape_space_dimensionality', inf), ...
%             struct('method', 'shortest path',   'weight_factor', 2, 'desired_shape_space_dimensionality', inf), ...
%             struct('method', 'regress',         'weight_factor', 0, 'desired_shape_space_dimensionality', inf), ...
%             struct('method', 'regress',         'weight_factor', 1, 'desired_shape_space_dimensionality', inf), ...
%             struct('method', 'regress',         'weight_factor', 2, 'desired_shape_space_dimensionality', inf), ...
%             struct('method',  'shortest path',   'weight_factor', 0, 'desired_shape_space_dimensionality', 7), ...
%             struct('method', 'shortest path',   'weight_factor', 1, 'desired_shape_space_dimensionality', 7), ...
%             struct('method', 'shortest path',   'weight_factor', 2, 'desired_shape_space_dimensionality', 7), ...
%             struct('method', 'regress',         'weight_factor', 0, 'desired_shape_space_dimensionality', 7), ...
%             struct('method', 'regress',         'weight_factor', 1, 'desired_shape_space_dimensionality', 7), ...
%             struct('method', 'regress',         'weight_factor', 2, 'desired_shape_space_dimensionality', 7), ...
%             };
%                 struct('method', 'nystrom-euclidean', 'desired_shape_space_dimensionality', 7), ...
%             struct('method', 'landmark-kernel-svd', 'desired_shape_space_dimensionality', 7), ...
%             struct('method', 'landmark-distance-svd', 'desired_shape_space_dimensionality', 7)

sampling = {'column'};

niter = 10;

if ~exist('distances_complete', 'var')
    %simulate some gaussian distrubuteddata
    
    ndims = 10;
    npts = 250;
    pts = mvnrnd(zeros(1, ndims), eye(10), npts);
    
    distances_complete = squareform(pdist(pts));
else
    distances_complete(eye(size(distances_complete))>0) = nan;
    
    keepinds = find(~all(isnan(distances_complete),2));
    
    distances_complete = distances_complete(keepinds,keepinds);
    distances_complete(eye(size(distances_complete))>0) = 0;
    
    nancounts = sum(isnan(distances_complete),2);
    
    while any(nancounts)
        [~, ind] = max(nancounts);
        
        distances_complete(ind,:) = [];
        distances_complete(:, ind) = [];
        keepinds(ind) = [];
        
        nancounts = sum(isnan(distances_complete),2);
    end
    
    npts = size(distances_complete,1);
end

eigval_best_temp = eig(distances_complete);
eigval_best_temp = eigval_best_temp / max(eigval_best_temp);
eigval_complete = zeros(1, npts);
eigval_complete(1:length(eigval_best_temp)) = eigval_best_temp(end:-1:1);

method = struct('method',  'regress',   'weight_factor', 0, 'desired_shape_space_dimensionality', inf);

posfile = [tempdir filesep 'posfile.mat'];
if ~exist(posfile, 'file')
    [positions_complete, ~, ~] = embed_partial_distance_matrix(distances_complete, method);
    save(posfile, 'positions_complete')
else
    load(posfile)
end

distances_euclidean = squareform(pdist(positions_complete));

for k = 1:niter
    
    
    for i = 1:length(methods)
        
        file_root = [tempdir filesep 'temp' num2str(k) '_' num2str(i) ];
        matfile = [file_root '.mat'];
        tmpfile = [file_root '.tmp'];
        
        if ~exist(matfile, 'file') && ~exist(tmpfile, 'file')
            system(['touch "' tmpfile '"'])
            
            disp(methods{i}.method);
            
            this_data = sparse(npts, npts);
            this_data(logical(eye(size(this_data)))) = true;
            
            rand_columns = randperm(npts);
            
            try
                
                for j = 1:npts
                    disp([num2str(j) filesep num2str(npts)]);
                    
                    inds = find(~all(this_data,1));
                    
                    if j == 1
                        endind = 1;
                    else
                        endind = 0;
                    end
                    
                    switch sampling{1}
                        case 'random' %this is currently not implemented correctly
                            npts = sum(~this_data(inds(1:endind),:) == 0);
                            
                            unsampled_points = find(this_data == 0);
                            
                            inds = randperm(length(unsampled_points));
                            
                            this_data(inds(npts)) = true;
                        case 'column'
                            rand_column = rand_columns(j);
                            
                            this_data(rand_column,:) = true;
                            this_data(:,rand_column) = true;
                        otherwise
                    end
                    
                    fract_data(j) = sum(this_data(:)>0)/ (npts^2);
                    
                    distances_temp = distances_complete;
                    distances_temp(this_data == 0) = nan;
                    
                    %                 methods{i} = ml_initparam(methods{i}, struct
                    
                    [positions{j}, ~, embedding_info{j}] = embed_partial_distance_matrix(distances_temp, methods{i});
                    
                end
                
                save(matfile, 'positions', 'fract_data', 'embedding_info', 'rand_columns');
                system(['rm ' tmpfile]);
            catch
                disp(['Method ' methods{i}.method ' failed.']);
            end
            
        end
    end
end

eignorm_vs_prev_all = ones(npts, length(methods), niter)* -1;
eignorm_vs_best_all = ones(npts, length(methods), niter)* -1;

matnorm_vs_prev_all = ones(npts, length(methods), niter)* -1;
matnorm_vs_best_all = ones(npts, length(methods), niter)* -1;

pos_vs_prev_all = ones(npts, length(methods), niter)* -1;
pos_vs_best_all = ones(npts, length(methods), niter)* -1;

corecoeff_vs_best_all = ones(npts, length(methods), niter)* -1;
corecoeff_euc_vs_best_all = ones(npts, length(methods), niter)* -1;

stress = ones(npts, length(methods), niter)* -1;
bic = ones(npts, length(methods), niter)* -1;

for k = 1:niter
    for i = 1:length(methods)
        matfile = [tempdir filesep 'temp' num2str(k) '_' num2str(i) '.mat'];
        if exist(matfile, 'file')
            try
                load(matfile)
            catch
                continue
            end
            
            for j = 1:npts
                
                dists = squareform(pdist(real(positions{j})));
                eigval_tmp = eig(dists);
                eigval_tmp = eigval_tmp / max(eigval_tmp);
                
                eigval = zeros(1, npts);
                eigval(1:length(eigval_tmp)) = eigval_tmp(end:-1:1);
                
                if ~exist('dists_prev', 'var')
                    dists_prev = dists;
                    eigval_prev = eigval;
                    positions_prev = positions{j};
                    
                end
                
                eignorm_vs_prev_all(j,i,k) = sum((eigval(:) - eigval_prev(:)).^2);
                eignorm_vs_best_all(j,i,k) = sum((eigval(:) - eigval_complete(:)).^2);
                
                matnorm_vs_best_all(j,i,k) = mean((dists(:) - distances_complete(:)).^2);
                
                pos_vs_prev_all(j,i,k) = procrustes_l2(positions{j}, positions_prev);
                pos_vs_best_all(j,i,k) = procrustes_l2(positions{j}, positions_complete);
                
                corr = corrcoef(dists(:), distances_complete(:));
                corecoeff_vs_best_all(j,i,k) = corr(1,2).^2;
                
                corr_euc = corrcoef(dists(:), distances_euclidean(:));
                corecoeff_euc_vs_best_all(j,i,k) = corr_euc(1,2).^2;
                
                stress(j,i,k) = embedding_info{j}.stress;
                
                upper_tri = triu(ones(npts,npts),1)>0;
                
                dists_vs_prev(j,i,k) = sum((dists(upper_tri) - dists_prev(upper_tri)).^2);
                
                legalmat = false(npts, npts);
                legalmat(rand_columns(1:j), :) = true;
                legalmat(:,rand_columns(1:j)) = true;
                legalmat = triu(legalmat,1);
                
                sample_var = var(distances_complete(legalmat));
                sum_error = sum((dists(legalmat) - distances_complete(legalmat)).^2);
                
                ndims = size(positions{j},2);
                bic(j,i,k) = 1/sample_var.^2*sum_error +  ndims*npts*log((npts*(npts-1))/2);
                
                dists_prev = dists;
                eigval_prev = eigval;
                positions_prev = positions{j};
            end
        end
    end
end



keepinds = ~all(all(matnorm_vs_best_all == -1,3),1);



matnorm_vs_best = matnorm_vs_best_all(:,keepinds,:);
eignorm_vs_best = eignorm_vs_best_all(:,keepinds,:);
eignorm_vs_prev = eignorm_vs_prev_all(:,keepinds,:);
matnorm_vs_prev = matnorm_vs_prev_all(:,keepinds,:);
pos_vs_best = pos_vs_best_all(:,keepinds,:);
pos_vs_prev = pos_vs_prev_all(:,keepinds,:);

corecoeff_vs_best = corecoeff_vs_best_all(:,keepinds,:);
corecoeff_euc_vs_best = corecoeff_euc_vs_best_all(:,keepinds,:);


stress = stress(:,keepinds,:);



xdat = (1:npts)./npts;

colormap(jet)

% matnorm_vs_best = cat(2, matnorm_vs_best(:,end-1:-1:1,:), matnorm_vs_best(:,end,:));
% matnorm_vs_best = matnorm_vs_best/numel(distances_complete);
%
% corecoeff_vs_best = cat(2, corecoeff_vs_best(:,end-1:-1:1,:), corecoeff_vs_best(:,end,:));
% corecoeff_euc_vs_best = cat(2, corecoeff_euc_vs_best(:,end-1:-1:1,:), corecoeff_euc_vs_best(:,end,:));

h(1) = figure('color', 'w');
plot_mean_and_stdev(xdat, matnorm_vs_best);
% title('MSE embedding distance matrix vs complete distance matrix', 'FontSize', 8)
xlabel('Fraction of observed data ');
ylabel('MSE vs complete data ', 'FontSize', 8);
axis square
axis tight

% h_legend = legend(h, functioning_methods);
% set(h_legend, 'FontSize', 6)


h(2) = figure('color', 'w');
plot_mean_and_stdev(xdat, 1-corecoeff_euc_vs_best);
% title('correlation coefficient vs complete data')
xlabel('Fraction of observed data ');
ylabel('1 - R^2');
% h_legend = legend(h, functioning_methods, 'Location', 'SouthWest');
% set(h_legend, 'FontSize', 6)

axis([0,1,0,1])
axis square
axis tight


for i = 1:length(h)
    
    set(get(h(i), 'CurrentAxes'), 'FontSize', 8)
    set(get(h(i), 'CurrentAxes'),'position',[0.15 0.15 .8 .8],'units','normalized')
    %     set(get(h(i), 'CurrentAxes'), 'Axes', 'tight')
    set(h(i), 'PaperUnits', 'Inches')
    set(h(i), 'PaperSize', [1.66, 1.66])
end
print([savedir filesep 'dmat_vs_complete.eps'], h(1), '-depsc2', '-r600')
print([savedir  filesep 'corrcoef_euc_vs_complete.eps'], h(2), '-depsc2', '-r600')





%
% [~, ind] = sort(sum(sum(matnorm_vs_best,3),1));
% functioning_methods(ind)'
%
%
% figure('color', 'w')
% h = plot_mean_and_stdev(xdat, pos_vs_best);
% title('L2 positions vs complete data')
% xlabel('Fraction of completed columns')
% ylabel('L2 vs complete data');
% legend(h, functioning_methods)
%
% if exist(savedir, 'dir')
%     saveas(gcf, [savedir filesep mfilename '_positions_vs_complete.tif'], 'tif')
% end
%
%
% figure('color', 'w')
% h = plot_mean_and_stdev(xdat, eignorm_vs_best);
% title('L2 eigenvalues vs complete data')
% xlabel('Fraction of completed columns')
% ylabel('L2 vs complete data');
% legend(h, functioning_methods)
%
% if exist(savedir, 'dir')
%     saveas(gcf, [savedir filesep mfilename '_eig_vs_complete.tif'], 'tif')
% end
%
% figure('color', 'w')
% h = plot_mean_and_stdev(xdat, matnorm_vs_prev);
% title('L2 distance matrix vs previous iteration')
% xlabel('Fraction of completed columns')
% ylabel('L2 vs complete data');
% legend(h, functioning_methods)
%
% if exist(savedir, 'dir')
%     saveas(gcf, [savedir filesep mfilename '_dmat_vs_prev.tif'], 'tif')
% end
%
% figure('color', 'w')
% h = plot_mean_and_stdev(xdat, eignorm_vs_prev);
% hold on
% title('L2 eigenvaluse vs previous iteration')
% xlabel('Fraction of completed columns')
% ylabel('L2 vs complete data');
% legend(h, functioning_methods)
%
% if exist(savedir, 'dir')
%     saveas(gcf, [savedir filesep mfilename '_eig_vs_prev.tif'], 'tif')
% end
%
%
% figure('color', 'w')
% h = plot_mean_and_stdev(xdat, pos_vs_prev);
% title('L2 positions vs previous iteration')
% xlabel('Fraction of completed columns')
% ylabel('L2 vs complete data');
% legend(h, functioning_methods)
%
% if exist(savedir, 'dir')
%     saveas(gcf, [savedir filesep mfilename '_pos_vs_prev.tif'], 'tif')
% end
%

end

function plot_mean_and_stdev(xdat, ydat, linewidth)
if ~exist('linewidth', 'var') | isempty(linewidth)
    linewidth = 1;
end

m = mean(ydat,3);
st = std(ydat,[], 3);

colors = jet(size(m,2))*0.8;

h = zeros(1, size(m,2));

for i = 1:size(m,2)
    
    plot(xdat, (m(:,i)), 'color', colors(i,:), 'LineWidth', linewidth);
    hold on
    plot(xdat, (m(:,i)+st(:,i)), '--', 'color', colors(i,:), 'LineWidth', linewidth)
    plot(xdat, (m(:,i)-st(:,i)), '--', 'color', colors(i,:), 'LineWidth', linewidth)
    
    
end
%     axis([0, 1, prctile(real(log10(m(m~=0))), 0.5), prctile(real(log10(m(m~=0))), 99.5)])
end