function [reconstructed_matrix] = reconstruct_distance_matrix2(known_distances, desired_dimensionality, options)
  % Directly solve for all entries in a distance matrix that is partially observed. In known_distances, entries are NaN when unknown. At least desired_dimensionality + 2 columns are assumed to be completely observed.
  %
  % Dependencies:
  %     Low-rank optimization for distance matrix completion
  %     Authors: B. Mishra, G. Meyer and R. Sepulchre
  %     http://www.montefiore.ulg.ac.be/~mishra/softwares/distCompletion.html
  %
  % Uses code C&P'd from the example
  %
  % By Gregory R. Johnson 8/19/2013 gj@andrew.cmu.edu

%   default_options = struct();
%   default_options.true_distance_matrix = [];
%   % default_options.use_all_known_columns = false;
%   default_options.use_all_known_columns = true;
%   default_options.use_squared_distances = true;
%   default_options.zero_negative_results = true;
%   default_options.debug = false;


    % Operate on the squared distance matrix:
    known_distances = known_distances.^2;
  
  n = size(known_distances,1);
  % Rank of the solution
    r = 10;
    % Starting approximation rank
    p = 1;

    % Parameters
    params.pmax = r;
    params.tol = 1e-3;
    params.vtol = 1e-3;
    params.verb = false; % "Bark" only when asked

  
    [I,J] = meshgrid(1:size(known_distances,1));

    distinds = intersect(find(triu(known_distances)), find(~isnan(known_distances)));

    methodFun = @tr_dist_completion;

    Y0 = randn(n, p);
    
    [Y infos] = lowrank_dist_completion(methodFun,I(distinds),J(distinds),known_distances(distinds),Y0,params);
  
    reconstructed_matrix = squareform(pdist(Y).^2);
    
    reconstructed_matrix = sqrt(reconstructed_matrix);
end
