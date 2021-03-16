% function [ output_args ] = proj_error_eig_test( model, savedir )


distances_incomplete = model.cellShapeModel.distances_incomplete;

distances_incomplete(logical(eye(size(distances_incomplete)))) = nan;

keepinds = find(~all(isnan(distances_incomplete),1));
    
distances_complete = distances_incomplete(keepinds, keepinds);
distances_complete(logical(eye(size(distances_complete)))) = 0;

landmarkinds = find(all(~isnan(distances_complete),1));
nlandmarks = size(landmarkinds,2);

distances_complete = distances_complete(landmarkinds, landmarkinds);


niter = 100;


figure, hold on

%for each iteration

embedding_err_prev              = cell(1, niter);

h = zeros(niter, nlandmarks);
p = zeros(niter, nlandmarks);

spectrum_err = zeros(niter, nlandmarks);

for i = 1:niter
    i
    landmarkpool = 1:nlandmarks;
    
    newind = randsample(nlandmarks, 1);
    testpool = landmarkpool(newind);
    landmarkpool(newind) = [];

    %starting at 2 landmarks
    for j = 2:nlandmarks
%       j
        %add a new landmark
        newind = randsample(nlandmarks-(j-1), 1);
        testpool(j) = landmarkpool(newind);
        landmarkpool(newind) = [];
        

        test_distances = distances_complete(testpool, testpool);
        
        tic

        [best_arrangement{i,j}, mass_matrix{i,j}, embedding_info{i,j}] = embed_partial_distance_matrix(test_distances, struct('method',  'nystrom-euclidean', 'force_positive_definiteness', false));
        
       
%         if j >= 3
%             ndims = size(best_arrangement{i,j},2);
%             n_dims_prev = size(best_arrangement{i,j-1},2);
%             
%             %exclude the most recently added point
%             best_arrangement_temp = best_arrangement{i,j}(1:end-1,:);
%             
%             spectrum = diag(embedding_info{i,j}.extended_kernel_eigenvalues);
%             spectrum = spectrum/sum(spectrum);
%             spectrum_prev = [diag(embedding_info{i,j-1}.extended_kernel_eigenvalues)];
%             spectrum_prev = spectrum_prev/sum(spectrum_prev);
%             
%             if ndims > n_dims_prev
%                 best_arrangement_prev = [best_arrangement_prev, zeros(j-1, ndims - n_dims_prev)];
%             elseif ndims < ndims_approx
%                 best_arrangement_temp = [best_arrangement_temp, zeros(j-1, n_dims_prev - ndims)];
%             end
%             
%             embedding_err_prev{i}(j) = procrustes(best_arrangement_temp, best_arrangement_prev);
%             
%             spectrum_err(i,j) = sum((spectrum(1:end-1) - spectrum_prev).^2)/ length(spectrum_prev);
%             [h(i,j), p(i,j)] = kstest2(spectrum, spectrum_prev);
%             
%             
%         end
%        
%         best_arrangement_prev = best_arrangement{i,j};
    end
end

for i = 1:size(embedding_info,1)
    for j = 3:size(embedding_info,2)
        spectrum = diag(embedding_info{i,j}.unfiltered_result_eigenvalues);
        spectrum = spectrum/sum(abs(spectrum));
        
        spectrum_prev = [diag(embedding_info{i,j-1}.unfiltered_result_eigenvalues)];
        spectrum_prev = spectrum_prev/sum(abs(spectrum_prev));
        
        [~, ind] = min(abs(spectrum_prev));
        
        if sign(spectrum_prev(ind)) == 1
            spectrum_prev = [spectrum_prev(1:ind); 0;spectrum_prev(ind+1:end)];
        else
            spectrum_prev = [spectrum_prev(1:ind-1); 0; spectrum_prev(ind:end)];
        end
        
        spectrum_err(i,j) = sum((spectrum - spectrum_prev).^2)/ length(spectrum_prev);
        [h(i,j), p(i,j)] = kstest2(spectrum, spectrum_prev);
        
    end
end

if exist('savedir','var')
    if ~exist(savedir, 'dir')
        mkdir(savedir)
    end
    
    save([savedir filesep 'proj_error_nystrom_simulate_vs_prev'])
end

embedding_err_prev_mean = mean(vertcat(embedding_err_prev{:}));
embedding_err_prev_mean = embedding_err_prev_mean(:,3:end);
embedding_err_prev_std = std(vertcat(embedding_err_prev{:}));
embedding_err_prev_std = embedding_err_prev_std(:,3:end);

figure('color', 'w')
hold on
plot(3:nlandmarks, embedding_err_prev_mean)
plot(3:nlandmarks, embedding_err_prev_mean+embedding_err_prev_std, '--')
plot(3:nlandmarks, embedding_err_prev_mean-embedding_err_prev_std, '--')
xlabel('# landmarks')
ylabel('SSE')
title('Embedding error vs previous')

% end

