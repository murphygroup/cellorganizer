function [ err ] = proj_error_eig( distances_incomplete )
%Procedure:

% Remove a landmark column
%     Compute an embedding using nystrom
%     compute normalizes SSE of eigenvalues


distances_incomplete(logical(eye(size(distances_incomplete)))) = nan;

keepinds = find(~all(isnan(distances_incomplete),1));
    
distances_complete = distances_incomplete(keepinds, keepinds);
distances_complete(logical(eye(size(distances_complete)))) = 0;

landmarkinds = find(all(~isnan(distances_complete),1));
nlandmarks = size(landmarkinds,2);


eigvals = zeros(nlandmarks, nlandmarks-1);

for i = 1:nlandmarks
    keepinds = true(1,length(landmarkinds));
    keepinds(i) = 0;
   
    [~, eigval, ~] = svd(distances_complete(:,landmarkinds(keepinds)));
    
    eigvals(i,:) = diag(eigval).^2;
    
end

sqerr = pdist(eigvals).^2;

err = sum(sqerr)/length(sqerr)/(nlandmarks-1);

end

