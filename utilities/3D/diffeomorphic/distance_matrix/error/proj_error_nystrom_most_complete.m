function [ output_args ] = proj_error_eig_test( model )


distances_incomplete = model.cellShapeModel.distances_incomplete;

distances_incomplete(logical(eye(size(distances_incomplete)))) = nan;

keepinds = find(~all(isnan(distances_incomplete),1));
    
distances_complete = distances_incomplete(keepinds, keepinds);
distances_complete(logical(eye(size(distances_complete)))) = 0;

landmarkinds = find(all(~isnan(distances_complete),1));
nlandmarks = size(landmarkinds,2);

distances_complete = distances_complete(landmarkinds, landmarkinds);


for i = 1:nlandmarks
    i
    combinations = combinator(nlandmarks, i, 'c');
    
    for j = 1:size(combinations,1)
        errs{i}(j) = proj_error_nystrom(distances_complete(combinations(j,:), combinations(j,:)));
    end
    
    err_mean(i) = mean(errs{i});
    err_std(i) = std(errs{i});

end

figure, hold on
plot(1:nlandmarks, err_mean, 'r');
plot(1:nlandmarks, err_mean+err_std, '--r');
plot(1:nlandmarks, err_mean-err_std, '--r');
xlabel('number landmarks')
ylabel('error')


end

