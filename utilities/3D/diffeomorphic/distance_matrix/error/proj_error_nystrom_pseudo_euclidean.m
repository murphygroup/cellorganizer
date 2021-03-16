function [ approximate_arrangement, approximate_arrangement_mass_matrix, hold_out ] = proj_error_nystrom( distances_incomplete )
%Procedure:

% Remove a landmark column
%     Compute an embedding using nystrom
%     compute error between the true distances and the embedded distances


    distances_incomplete(logical(eye(size(distances_incomplete)))) = nan;

    keepinds = find(~all(isnan(distances_incomplete),1));

    distances_complete = distances_incomplete(keepinds, keepinds);
    distances_complete(logical(eye(size(distances_complete)))) = 0;

    landmarkinds = find(all(~isnan(distances_complete),1));
    nlandmarks = size(landmarkinds,2);


    distances_complete = distances_complete(landmarkinds, landmarkinds);

%     err_all = zeros(1, nlandmarks);
% 
%     [best_arrangement, ~, ~] = embed_partial_distance_matrix(distances_complete, struct('method',  'nystrom-euclidean', 'force_positive_definiteness', true));
    
%     ndims = size(best_arrangement,2);
%     
%     [~,b,~] = eig(distances_complete);
%     eigval = diag(b);
%     
%     npts = size(distances_complete,1);
%     
    approximate_arrangement = cell(1, nlandmarks);
    approximate_arrangement_mass_matrix = cell(1, nlandmarks);
    holdout = zeros(1, nlandmarks);
    
    for i = 1:nlandmarks
%         best_arrangement_temp = best_arrangement;
%         
        dists = distances_complete;
        dists(:,i) = 0;

        [approximate_arrangement{i}, approximate_arrangement_mass_matrix{i}, approximate_arrangement_embedding_info] = embed_partial_distance_matrix(dists, struct('method',  'nystrom-euclidean', 'force_positive_definiteness', true));
        hold_out(i) = i;
%         [~,b,~] = eig(approximate_arrangement_distances);
%         
%         eigval_hoo = diag(b);
%         
%         err_all(i) = sum((eigval_hoo - eigval).^2)/(nlandmarks-1);

    end

end



