function [ output_args ] = proj_error_eig_test( model, savedir )


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

embedding_err_complete_iter     = cell(niter, nlandmarks);
embedding_err_per_point_iter    = cell(niter, nlandmarks);
distances_err_complete_iter     = cell(niter, nlandmarks);
distances_err_per_point_iter    = cell(niter, nlandmarks);

embedding_err_prev              = cell(1, niter);


savefile = [savedir filesep 'proj_error_nystrom_simulate.mat'];

if ~exist(savefile, 'file')

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
            %get the hold-out-one embedddings

            [best_arrangement{i,j}, ~, ~] = embed_partial_distance_matrix(test_distances, struct('method',  'nystrom-euclidean', 'force_positive_definiteness', true));
            [ approximate_arrangement, approximate_arrangement_mass_matrix, hold_out] = proj_error_nystrom_pseudo_euclidean(test_distances);

            npts = size(best_arrangement{i,j},1); % equalivalently npts = j
            ndims = size(best_arrangement{i,j},2);

            embedding_err_complete = zeros(1, length(hold_out));
            embeddingt_err_per_point = zeros(1, length(hold_out));
            distances_err_complete = zeros(1, length(hold_out));
            distances_err_per_point = zeros(1, length(hold_out));



            for k =  1:length(hold_out)

                approximate_distance_function = @(given_vectors1, given_vectors2)mass_distance_function(given_vectors1, given_vectors2, approximate_arrangement_mass_matrix{j});

                approximate_arrangement_distances = squareform(pdist(approximate_arrangement{j}, approximate_distance_function));

                ndims_approx = size(approximate_arrangement{k},2);

                best_arrangement_temp = best_arrangement{i,j};

                if ndims > ndims_approx
                    approximate_arrangement{k} = [approximate_arrangement{k}, zeros(npts, ndims - ndims_approx)];
                elseif ndims < ndims_approx
                    best_arrangement_temp = [best_arrangement_temp, zeros(npts, ndims_approx - ndims)];
                end

                [embedding_err_complete(k), approximate_arrangement_aligned] = procrustes(best_arrangement_temp,approximate_arrangement{k});

                embeddingt_err_per_point(k) = sum((approximate_arrangement_aligned(hold_out(k),:) - best_arrangement_temp(hold_out(k),:)).^2);

                distances_err_complete(k) = sum((approximate_arrangement_distances(:) - test_distances(:)).^2);

                distances_err_per_point(k) = sum((approximate_arrangement_distances(hold_out(k),:) - test_distances(hold_out(k),:)).^2);

            end


            if j >= 3
                n_dims_prev = size(best_arrangement_prev,2);

                %exclude the most recently added point
                best_arrangement_temp = best_arrangement{i,j}(1:end-1,:);

                if ndims > n_dims_prev
                    best_arrangement_prev = [best_arrangement_prev, zeros(npts-1, ndims - n_dims_prev)];
                elseif ndims < ndims_approx
                    best_arrangement_temp = [best_arrangement_temp, zeros(npts-1, n_dims_prev - ndims)];
                end

                embedding_err_prev{i}(j) = procrustes(best_arrangement_temp, best_arrangement_prev);
            end

            embedding_err_complete_iter{i,j} = embedding_err_complete;
            embedding_err_per_point_iter{i,j} = embeddingt_err_per_point;
            distances_err_complete_iter{i,j} = distances_err_complete;
            distances_err_per_point_iter{i,j} = distances_err_per_point;

            best_arrangement_prev = best_arrangement{i,j};

        end

    end

    % plot(1:nlandmarks, err_mean, 'r');
    % plot(1:nlandmarks, err_mean+err_std, '--r');
    % plot(1:nlandmarks, err_mean-err_std, '--r');
    % xlabel('number landmarks')
    % ylabel('error')
    if exist('savedir','var')
        if ~exist(savedir, 'dir')
            mkdir(savedir)
        end

        save(savefile)
    end

else
    load(savefile)
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

end

