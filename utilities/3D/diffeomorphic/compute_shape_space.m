function [result, success] = compute_shape_space(...
  distances, param, save_filename)

warning('This function is depreicated. Please use ''embed_partial_distance_matrix.m''')

%grj 9/22/13 - force to minimum of 2 dimensions
%grj 11/26/14 - removed the convex hull calculation

  % [result, success] = compute_shape_space(...
  %   distances, numdims, explained_variance_threshold, save_filename)
  %
  % Given a distance matrix of a set of points, computes positions for those points using multidimensional scaling (y2) (the space in which those points live is called a shape space), its convex hull (convex_hull), and the Delaunay tessellation (tes) of those points.
  %
  % If explained_variance_threshold is an array with two elements, the second is the maximum allowable dimensionality.
  
  %error('Implementation yet unfinished below this line!')
  
  if ~exist('param', 'var') | isempty(param)
      param = struct();
  end
  
  param = ml_initparam(param, struct( ...
                    'maxdims', min([7,size(distances,1)]), ...
                    'explained_variance_threshold', 0.90, ...
                    'weight_factor', 0, ...
                    'just_positions', false ...
                    ));
  
  result = struct();
  result.positions = [];
  result.convex_hull = [];
  result.numdims = [];
  
  positions = [];     
  convex_hull = [];
  tes = [];    
  
%   success = false;
%   % Check if work has been done
%   [can_start, final_name, final_exists, tmpfile] = chunk_start(save_filename);
%   
% %   if (~can_start)
% %     stack = dbstack();
% %     % warning('Skipping "%s" at line %d in file %s.', save_filename, stack(1).line, stack(1).file)
% %     fprintf('compute_shape_space: Skipping "%s" at line %d in file %s.\n', save_filename, stack(1).line, stack(1).file)
% %     if final_exists
% %       load([save_filename '.mat'])
% %       success = true;
% %     else
% %       %error('Cannot continue, serial 1 being computed elsewhere.')
% %       error('  Being computed elsewhere.')
% %     end
%   if ~can_start && ~final_exists
%     error('  Being computed elsewhere.')
%   else
    
    % keyboard
    
    
    if param.weight_factor == 0
        [y, e] = cmdscale( distances );
    else
        [y] = mdscale(distances, size(distances,1), 'Weights', 1./(distances.^param.weight_factor), 'Options', statset('Display', 'iter'));
        e = eig(y*y');
        e = e(end:-1:1);
    end
    
    % e

    % Determine shape space dimensionality using positive eigenvalues
    % (too bad some of the negative ones can be large):
%     explained_variances = [];
%     if param.maxdims == 0
      explained_variances = cumsum(e .* (e >= 0));
      explained_variances = explained_variances ./ explained_variances(end);
      sufficient_dimensionalities = ...
          find(explained_variances >= ...
               param.explained_variance_threshold); 
      numdims = sufficient_dimensionalities(1);
      numdims = min(numdims, param.maxdims);
     
%     end
    
    if numdims > size(y, 2)
      warning(['numdims = ', num2str(numdims), ' > size(y, 2) = ', num2str(size(y, 2))])
      numdims = size(y, 2);
    end
    
    if numdims == 1
        numdims = 2;
    end
    
    convex_hull = [];
    tes = [];
    if size(y,2) == 1
        positions = y;
%         convex_hull = nan;
        tes = nan;
    else 
        positions = y(:,1:numdims);

        if ~param.just_positions
            convex_hull = convhulln( positions );
            %crashes:
            %convex_hull = convhulln( y2, {'Qt','Qx', 'T1'} );
            tes = delaunayn(positions);
        end
    end
    
    if exist('save_filename', 'var') & ~isempty(save_filename)
        save([save_filename '.mat'], 'positions', 'convex_hull', 'tes', 'numdims', 'explained_variances');
    end
    %%%%%%%%%%%%%%%%
% 
%     % Save work
%     chunk_finish(tmpfile); 
%   end

  result.positions = positions;
  result.convex_hull = convex_hull;
  result.tessellation = tes;
  result.numdims = numdims;
  result.explained_variances = explained_variances;
  
  success = true;
