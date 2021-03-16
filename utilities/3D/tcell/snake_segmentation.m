function [result] = snake_segmentation(image, seed_points, seed_radii, options)
% Run Snake3D from Dirk de Kroon with spherical seeds.
%
% seed_points is Nx3 with XYZ ordering or a patch-style structure.
%
% If seed_points is a set of points, a spherical surface is placed at each seed point along with the corresponding radius in seed_radii, otherwise the given mesh is used. Snake segmentation then runs.
%
% 2012-09-01 tebuck: created.
% 2012-11-18 tebuck: added options.replicate_padding to prevent spurious edges at the image boundary.

  % input_options = options
  
  % Process options:
  default_options = struct(...
    'relative_convergence_criterion', 1e-2 ...
    , 'seed_subdivisions', 2 ...
    , 'replicate_padding', false ...
    ... % These are the defaults from Dirk de Kroon's Snake3D:
    , 'Verbose', false ...
    , 'Wline', 0.04 ...
    ,'Wedge', 2 ...
    ,'Sigma1', 2 ...
    ,'Sigma2', 2 ...
    ,'Alpha', 0.2 ...
    ,'Beta', 0.2 ...
    ,'Delta', 0.1 ...
    ,'Gamma', 1 ...
    ,'Kappa', 2 ...
    ,'Iterations', 100 ...
    ,'GIterations', 0 ...
    ,'Mu', 0.2 ...
    ,'Sigma3', 1 ...
    ,'Lambda', 0.8 ...
    );

  if ~exist('options', 'var')
    options = default_options; 
  else
    option_names = fieldnames(default_options); 
    for index = 1:length(option_names)
      current_option = option_names{index}; 
      if ~isfield(options, current_option)
        options = setfield(options, current_option, getfield(default_options, current_option)); 
      end
    end
  end
  %options
  
  all_seeds_mesh = [];
  if isstruct(seed_points)
    all_seeds_mesh = seed_points;
  else
    number_seed_points = size(seed_points, 1);
    if length(seed_radii) == 1
      seed_radii = ones(number_seed_points, 1) .* seed_radii; 
    end
    
    
    
    % Approximate unit sphere:
    % fprintf('>> CREATING SPHERE MESH\n')
    [v, f] = compute_semiregular_sphere(options.seed_subdivisions);
    seed_mesh = struct('vertices', v', 'faces', f');
    number_sphere_points = size(seed_mesh.vertices, 1);

    all_seeds_mesh = struct('vertices', zeros(0, 3), 'faces', zeros(0, 3)); 
    for seed_index = 1:number_seed_points
      current_seed_mesh = seed_mesh; 
      % size(current_seed_mesh.vertices .* seed_radii(seed_index))
      % size(repmat(seed_points(seed_index, :), number_sphere_points, 1))
      current_seed_mesh.vertices = ...
        current_seed_mesh.vertices .* seed_radii(seed_index) + ...
        repmat(seed_points(seed_index, :), number_sphere_points, 1);
      current_seed_mesh.faces = current_seed_mesh.faces + number_sphere_points * (seed_index - 1);
      all_seeds_mesh.vertices = [all_seeds_mesh.vertices; current_seed_mesh.vertices];
      all_seeds_mesh.faces = [all_seeds_mesh.faces; current_seed_mesh.faces];
      % all_seeds_mesh
    end
  end
  
  % Snake3D appears to use YXZ ordering:
  all_seeds_mesh.vertices = all_seeds_mesh.vertices(:, [2, 1, 3]);
  
  
  if options.replicate_padding
    % Snake3D uses symmetric padding, counteract this by replicate padding the image:
    % pad_size = ceil(3 * max([options.Sigma1, options.Sigma2, options.Sigma3]));
    % Should probably include GVF's GIterations here. 
    pad_size = ceil(3 * max([options.Sigma1, options.Sigma2, sqrt(options.Sigma3^2 * options.GIterations)]));
    image = padarray(image, pad_size * ones(1, 3), 'replicate');
    all_seeds_mesh.vertices = all_seeds_mesh.vertices + pad_size;
  end


  % clf
  % patch(all_seeds_mesh, 'facecolor', [1, .4, .1]); camlight; lighting phong; axis equal;
  
  % all_seeds_mesh
  % 'size(image)'
  % size(image)
  warning('off', 'snake:unknownoption');
  segmentation_mesh = Snake3D(image, all_seeds_mesh, options);
  
  
  if options.replicate_padding
    segmentation_mesh.vertices = segmentation_mesh.vertices - pad_size;
  end


  % Return to XYZ ordering:
  segmentation_mesh.vertices = segmentation_mesh.vertices(:, [2, 1, 3]);
  result = segmentation_mesh;
  
  % clf
  % patch(segmentation_mesh, 'facecolor', [1, .4, .1]); camlight; lighting phong; axis equal;
  
  