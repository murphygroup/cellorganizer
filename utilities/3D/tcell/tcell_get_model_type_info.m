function [t_cell_info] = tcell_get_model_type_info(t_cell_info, options)
  % 
  %
  % 2016-02-25 xruan: Copied from master_script_get_model_type_info.m.
  % % Dependencies:
  % % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
  % % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh

  
  default_options = struct();
  % % Allow an automatic parameter search (only continuous parameters for now):
  % default_options.model_parameters = struct();
  % default_options.model_parameters.standardized_voxels = struct();
  % default_options.model_parameterse.cylindrical_zernike = struct();
  % default_options.model_parameters.cylindrical_zernike.highest_order = 6;
  
    
  if ~exist('options', 'var')
    options = default_options; 
  else
    % options = process_options_structure(default_options, options);
    options = process_options_structure(default_options, options, [], true);
  end
  % options
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Copy variables to local namespace:
  
  
  master_script_options = t_cell_info.options;
  
  % illustration_colormap_to_use = t_cell_info.plotting_info.illustration_colormap_to_use;
  
  template_options = t_cell_info.template_info.template_options;
  template_image = t_cell_info.template_info.template_image;
  template_cell_radius = t_cell_info.template_info.template_cell_radius;
  template_synapse = t_cell_info.template_info.template_synapse;
  template_crop_function = t_cell_info.template_info.template_crop_function;
  
  template_occupied_x_slices = t_cell_info.template_info.template_occupied_x_slices;
  template_occupied_x_slice_pad_function = t_cell_info.template_info.template_occupied_x_slice_pad_function;
  template_occupied_y_slices = t_cell_info.template_info.template_occupied_y_slices;
  template_occupied_y_slice_pad_function = t_cell_info.template_info.template_occupied_y_slice_pad_function;
  approximate_synapse_slice = t_cell_info.template_info.approximate_synapse_slice;
  approximate_synapse_slice_template_image = t_cell_info.template_info.approximate_synapse_slice_template_image;
  approximate_cytoplasm_slices = t_cell_info.template_info.approximate_cytoplasm_slices;
  approximate_cytoplasm_slice_template_image = t_cell_info.template_info.approximate_cytoplasm_slice_template_image;
  template_x = t_cell_info.template_info.template_x;
  template_y = t_cell_info.template_info.template_y;
  template_z = t_cell_info.template_info.template_z;
  template_cylindrical_r = t_cell_info.template_info.template_cylindrical_r;
  template_cylindrical_theta = t_cell_info.template_info.template_cylindrical_theta;
  
  model_type_to_include = master_script_options.model_type_to_include;
  
  % convert the cellorganizer model type to tcell model type
  for i = 1 : numel(model_type_to_include)
      if strcmp('standardized_map_half-ellipsoid', model_type_to_include{i})
          model_type_to_include{i} = 'standardized_voxels';
      end
  end          
  
  % Sometimes causes errors with at least truncated_normal:
  % minimum_density = 0;
  minimum_density = 1e-3;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Create model types and associated functions:
  
  % xruan 02/25/2016
  model_types_all = {};
  % All the voxels inside the template:
  model_types_all{end + 1} = 'standardized_voxels';
  if false
    warning('>>>> HACK, only using one model type for debugging!')
  else
    % Mean intensity along the axis:
    model_types_all{end + 1} = 'marginal_x';
    % Intensity at the synapse:
    model_types_all{end + 1} = 'synapse_yz';
    % Horizontal marginal at the synapse:
    model_types_all{end + 1} = 'syn_mar_y';
    % Intensity at several slices near the synapse that contain much of the cytoplasm:
    model_types_all{end + 1} = 'fwd_cyto';
    % Two-parameter (peakedness and radius) radial intensity distribution in same region as fwd_cyto:
    % model_types{end + 1} = 'fwd_cyto_radial';
    % error('This isn''t ready yet')
    % radial_distribution_types = {'lognormal', 'beta', 'truncated_normal', 'blurred_uniform', 'constrained_beta'}
    radial_distribution_types = {'lognormal', 'beta', 'truncated_normal', 'constrained_beta', 'raw'};
    % radial_distribution_types = {'raw'}
    number_radial_distribution_types = length(radial_distribution_types);
    radial_distribution_model_types = strcat('fwd_cyto_radial_', radial_distribution_types);
    for radial_distribution_type_index = 1:number_radial_distribution_types
      model_types_all{end + 1} = radial_distribution_model_types{radial_distribution_type_index};
    end
  end
  
  for i = 1 : numel(model_type_to_include)
      if ~any(cellfun(@(x) strcmp(model_type_to_include{i}, x), model_types_all))
          error(sprintf('Model type %s does not exist', model_type_to_include{i}))
      end
  end
  
  model_types = model_type_to_include;
  number_model_types = length(model_types);
  model_representation_functions = cell(number_model_types, 1);
  model_reconstruction_functions = cell(number_model_types, 1);
  model_reduction_reconstruction_functions = cell(number_model_types, 1);
  model_reduction_representation_functions = cell(number_model_types, 1);
  model_original_representation_functions = cell(number_model_types, 1);
  model_original_reconstruction_functions = cell(number_model_types, 1);
  % model_should_compute_covariance = false(number_model_types, 1);
  % model_should_compute_chi_square = false(number_model_types, 1);
  % model_dimensionality_reduction_method_types = {'none', 'pca'};
  % If 'pca', compute PCA of all model vectors:
  model_dimensionality_reduction_methods = repmat({'none'}, [number_model_types, 1]);
  % model_dimensionality_reduction_method_parameters = cell(number_model_types, 1);
  model_dimensionality_reduction_method_parameters = repmat({struct()}, number_model_types, 1);
  model_dimensionality_reduction_method_transformations = cell(number_model_types, 1);
  model_illustration_view = repmat({'top'}, [number_model_types, 1]);
  
  blank_template_with_nans = nan(size(template_image));
  blank_template_with_nans(template_image > 0) = 0;
  blank_template_with_nans_inverse = nan(size(template_image));
  blank_template_with_nans_inverse(template_image == 0) = 0;
  
  
  % Useful values for fwd_cyto_radial_raw model:
  radial_raw_number_bins = round(template_cell_radius);
  % radial_raw_bin_centers = ((1:radial_raw_number_bins) - .5) / radial_raw_number_bins;
  % radial_raw_bin_edges = [-inf, (1:radial_raw_number_bins - 1) / radial_raw_number_bins, inf];
  % radial_raw_bin_edges = [-1, (1:radial_raw_number_bins - 1) / radial_raw_number_bins, 2];
  radial_raw_bin_edges = (0:radial_raw_number_bins) / radial_raw_number_bins;
  
  
  if false
    % Debug info:
    % % This works:
    a = randn(size(template_occupied_x_slices)); b = a(template_occupied_x_slices); b = template_occupied_x_slice_pad_function(b); whos a b, a(template_occupied_x_slices) - b(template_occupied_x_slices)
    % keyboard
  end
  
  function [result_density] = robust_blurred_uniform_pdf(gx, ga, gr, gs)
    % Produce a blurred uniform pdf that smoothly degrades into a normal when the parameters make it close to such:
    interpolation_factor = max(min((gs ./ gr - 9e1) ./ 1e1, 1), 0);
    % % Hack: interpolation_factor always > 0 so large values to erf calls below don't return 0 density:
    % interpolation_factor = max(interpolation_factor, 1e-2);
    if interpolation_factor < 1
      result_density = (1 ./ (2 .* (gr)) .* (erf(((ga + gr) - gx) ./ (gs .* sqrt(2))) - erf((ga - gx) ./ (gs .* sqrt(2)))));
      % result_density = (1 ./ (2 .* (gr)) .* (erf_(((ga + gr) - gx) ./ (gs .* sqrt(2))) - erf_((ga - gx) ./ (gs .* sqrt(2)))));
    end
    result_density = result_density .* (1 - interpolation_factor);
    if interpolation_factor > 0
      result_density = result_density + normpdf(gx, ga + gr ./ 2, gs) .* interpolation_factor;
    end
    % keyboard
  end
  
  function [result_parameters] = fwd_cyto_radial_model_representation(given_image, current_radial_distribution_type)
    valid_region = approximate_cytoplasm_slice_template_image;
    % radial_distribution_sampling = 'all';
    radial_distribution_sampling = 'even';
    switch radial_distribution_sampling
      case 'all'
        % Use all voxels:
        distribution_support = valid_region;
        if any(strcmp(radial_distribution_type, {'lognormal'}))
          distribution_support = distribution_support & (template_cylindrical_r > 0);
        end
        current_values = template_cylindrical_r(distribution_support);
        current_weights = given_image(distribution_support);
        % Weights should be based on intensity and angular size of each sample, not just intensity...
        current_weights = current_weights .* (1 ./ template_cylindrical_r(distribution_support));
      case 'even'
        % Sample uniformly:
        radial_sampling_rate = 1;
        % radial_sampling_rate = 2;
        % radial_sampling_rate = 4;
        angular_sampling_rate = 1;
        % angular_sampling_rate = 2;
        x_sampling_rate = 1;
        % x_sampling_rate = 2;
        given_image_nanned = given_image;
        given_image_nanned(~template_image) = nan;
        given_image_extrapolated = extrapolate_nans(given_image_nanned);
        r_limits = [0, 1];
        theta_limits = [0, 2 * pi];
        x_limits = [min(approximate_cytoplasm_slices), max(approximate_cytoplasm_slices)];
        r_number_values = ceil(template_cell_radius * radial_sampling_rate);
        theta_number_values = ceil(template_cell_radius * angular_sampling_rate * 2 * pi);
        x_number_values = ceil(length(approximate_cytoplasm_slices) * x_sampling_rate);
        % r_values = contrast_stretch(0:ceil(template_cell_radius * radial_sampling_rate));
        % if any(strcmp(radial_distribution_type, {'lognormal', 'beta', 'truncated_normal', 'constrained_beta'}))
        if true
          % For some reason fitdist doesn't let beta distros learn from values at 0 and 1:
          % r_values(1) = eps;
          % r_values(end) = 1 - eps;
          r_limits = [eps, 1 - eps];
        end
        % theta_values = (0:ceil(template_cell_radius * angular_sampling_rate * 2 * pi) - 1) ./ ceil(template_cell_radius * angular_sampling_rate * 2 * pi) * 2 * pi;
        % x_values = (0:ceil(length(approximate_cytoplasm_slices) * x_sampling_rate) - 1) ./ ceil(length(approximate_cytoplasm_slices) * x_sampling_rate) * length(approximate_cytoplasm_slices);
        
        % % Sample using a grid:
        % [r, theta, x] = ndgrid(r_values, theta_values, x_values);
        
        % Use Sobol sequence:
        sample_generator = sobolset(3);
        sampled_points = sample_generator(1:r_number_values * theta_number_values * x_number_values, :);
        % r = sampled_points(:, 1) .* range(r_values) + min(r_values);
        % theta = sampled_points(:, 2) .* (2 * pi);
        % x = sampled_points(:, 3) .* (2 * pi);
        r = sampled_points(:, 1) .* diff(r_limits) + min(r_limits);
        r = min(max(r, r_limits(1)), r_limits(2));
        theta = sampled_points(:, 2) .* diff(theta_limits) + min(theta_limits);
        theta = min(max(theta, theta_limits(1)), theta_limits(2));
        x = sampled_points(:, 3) .* diff(x_limits) + min(x_limits);
        x = min(max(x, x_limits(1)), x_limits(2));
        
        samples_x = x;
        samples_y = cos(theta) .* r * template_options.yr + template_options.yc;
        samples_z = sin(theta) .* r * template_options.zr + template_options.zc;
        % [samples_x, samples_y, samples_z] = sph2cart();
        current_values = r(:);
        current_weights = interp3(given_image_extrapolated, samples_y(:), samples_x(:), samples_z(:), 'linear');
        if false && master_script_options.debug_model_reconstruction
          % Debug info:
          % imagesc(reshape_contrast(given_image_extrapolated, -1)), axis image
          imagesc(squeeze(mean(given_image_extrapolated, 1))), axis image
          % colormap(pmkmp(1024, 'CubicL'))
          colorbar
          colormap(illustration_colormap_to_use)
          hold on
          scatter(samples_z(:), samples_y(:), [], samples_x(:))
          hold off
          beep, pause
        end
    end
    switch current_radial_distribution_type
      case 'normal_curve'
        curve_fit = fit(current_values, current_weights, 'gauss1');
        % http://www.mathworks.com/help/stats/examples/fitting-custom-univariate-distributions.html
        % Normalize:
        curve_fit.a1 = curve_fit.a1 / integrate(curve_fit, 1, 0);
        result_parameters = [curve_fit.a1, curve_fit.b1, curve_fit.c1];
        % error
        
            
      case {'truncated_normal', 'blurred_uniform', 'constrained_beta'}
        custom_pdf = [];
        custom_start = [];
        custom_lower = [];
        custom_upper = [];
        switch current_radial_distribution_type
          case 'truncated_normal'
            % Know range is [0, 1], so not free params:
            custom_pdf = @(gx, gm, gs)(normpdf(gx, gm, gs) .* ((gx >= 0) & (gx <= 1)) ./ diff(normcdf([0, 1], gm, gs))) .* (1 - minimum_density) + minimum_density;
            custom_start = [0.5, 1];
            custom_lower = [0, 1e-3];
            custom_upper = [1, 1e3];
          case 'constrained_beta'
            % Know range is [0, 1], so not free params:
            custom_pdf = @(gx, ga, gb)(betapdf(gx, ga, gb)) .* (1 - minimum_density) + minimum_density;
            custom_start = [2, 2];
            custom_lower = [1, 1];
            custom_upper = [1e3, 1e3];
          case 'blurred_uniform'
            custom_pdf = @(gx, ga, gr, gs)(robust_blurred_uniform_pdf(gx, ga, gr, gs)) .* (1 - minimum_density) + minimum_density;
            error('Doesn''t work well b/c giving zero density in places')
            if false
              % Debug info:
              plot((0:1000) ./ 1000, custom_pdf((0:1000) ./ 1000, .25, .5, 1e-1))
              beep
              pause
            end
            % error
            custom_start = [.25, .5, 5e-1];
            % custom_lower = [-inf, 0, 0];
            % custom_upper = [inf, inf, inf];
            custom_lower = [-2, 1e-2, 1e-3];
            custom_upper = [2, 4, 1e3];
          otherwise
            error('Unrecognized current_radial_distribution_type "%s"', current_radial_distribution_type)
        end
        
        current_weights = round(current_weights ./ double(max(current_weights)) * 1000);
        result_parameters = mle(current_values, 'frequency', current_weights, 'pdf', custom_pdf, 'start', custom_start, 'lower', custom_lower, 'upper', custom_upper);
        
        if false
          % Debug info:
          % result_parameters_cell_without_scale = num2cell(result_parameters(1:end - 1));
          result_parameters_cell_without_scale = num2cell(result_parameters);
          plot((0:1000) ./ 1000, custom_pdf((0:1000) ./ 1000, result_parameters_cell_without_scale{:}))
          beep
          % pause
          dbstack, keyboard
        end
        
        % error
        
      case 'raw'
        % Weighted histogram:
        [bin_counts, bin_assignments] = histc(current_values, radial_raw_bin_edges);
        result_parameters = accumarray(bin_assignments, current_weights);
        result_parameters = reshape(result_parameters, 1, []);
        result_parameters = result_parameters ./ sum(result_parameters);
        % result_parameters = result_parameters .* diff(radial_raw_bin_edges);
        result_parameters = result_parameters ./ diff(min(max(radial_raw_bin_edges, 0), 1));
        if false
          % Debug info:
          % warning('Finish implementing raw!!!')
          if length(result_parameters) == length(radial_raw_bin_edges)
            whos radial_raw_bin_edges result_parameters
            beep
            dbstack
            % pause
            keyboard
          end
        end

      otherwise
        % Frequency needs to be integer-valued (though lognfit doesn't require that, odd...):
        current_weights = round(current_weights ./ double(max(current_weights)) * 1000);
        distribution_fit = fitdist(current_values, current_radial_distribution_type, 'frequency', current_weights);
        % Normalize so probability mass in valid region is 1:
        distribution_cdf_in_range = cdf(distribution_fit, 1) - cdf(distribution_fit, 0);
        result_parameters = getfield(distribution_fit, 'Params');
    end
    
    result_parameters = [result_parameters, mean(given_image(valid_region))];
    
    if ~false && master_script_options.debug_model_reconstruction
      % Debug info:
      % [current_intervals, current_histogram] = histwc(current_values, current_weights, 20);
      % bar(current_histogram, current_intervals)
      % % scatter(current_values, current_weights)
      % beep, pause
      % KDE for pdf:
      kde_fit = fitdist(current_values, 'kernel', 'frequency', current_weights, 'width', .025);
      % kde_fit = fitdist(current_values, 'kernel', 'frequency', current_weights, 'support', [0, 1], 'width', .05);
      switch current_radial_distribution_type
        case 'normal_curve'
          % plot(curve_fit)
          % error
          plot((0:1000) ./ 1000, [curve_fit((0:1000) ./ 1000).'; kde_fit.pdf((0:1000) ./ 1000)]), legend({current_radial_distribution_type, 'kernel'})
        case {'truncated_normal', 'blurred_uniform', 'constrained_beta'}
          result_parameters_cell_without_scale = num2cell(result_parameters(1:end - 1));
          plot((0:1000) ./ 1000, [custom_pdf((0:1000) ./ 1000, result_parameters_cell_without_scale{:}); kde_fit.pdf((0:1000) ./ 1000)]), legend({current_radial_distribution_type, 'kernel'})
        case 'raw'
          % bar(radial_raw_bin_edges, [result_parameters(1:end - 1), 0.], 'histc')
          bar(radial_raw_bin_edges(1:end - 1), result_parameters(1:end - 1), 'histc')
        otherwise
          % 22:53 2014-06-08 beta fit looking horrible given num of samples... should be uniform :(
        plot((0:1000) ./ 1000, [distribution_fit.pdf((0:1000) ./ 1000); kde_fit.pdf((0:1000) ./ 1000)]), legend({current_radial_distribution_type, 'kernel'})
      end
      title('Radial intensity probability distribution', 'Interpreter', 'none')
      beep
      pause
      % keyboard
    end
    
    result_parameters = reshape(result_parameters, 1, []);
  end
  
  
  function [result_densities] = fwd_cyto_radial_model_density(given_parameters, current_radial_distribution_type, given_radial_coordinates)
    % Untested portion of fwd_cyto_radial_model_reconstruction:
    current_valid_region = (given_radial_coordinates >= 0) & (given_radial_coordinates <= 1);
    current_parameters_cell_without_scale = num2cell(given_parameters(1:end - 1));
    % Hack values so densities aren't actually computed at 0 or 1:
    given_radial_coordinates = min(given_radial_coordinates, 1 - eps * 1e10);
    given_radial_coordinates = max(given_radial_coordinates, eps * 1e10);
    switch current_radial_distribution_type
      case 'normal_curve'
        curve_fit = cfit(fittype('gauss1'), current_parameters_cell_without_scale{:})
        result_densities = curve_fit(given_radial_coordinates(current_valid_region)).';
        
      case {'truncated_normal', 'blurred_uniform', 'constrained_beta', 'raw'}
        % Don't duplicate code like this!
        switch current_radial_distribution_type
          case 'truncated_normal'
            % Know want truncated to range [0, 1], so not free params:
            % custom_pdf = @(gx, gm, gs)normpdf(gx, gm, gs) .* ((gx >= 0) & (gx <= 1)) ./ diff(normcdf([0, 1], gm, gs));
            custom_pdf = @(gx, gm, gs)(normpdf(gx, gm, gs) .* ((gx >= 0) & (gx <= 1)) ./ diff(normcdf([0, 1], gm, gs))) .* (1 - minimum_density) + minimum_density;
          case 'blurred_uniform'
            % custom_pdf = @robust_blurred_uniform_pdf;
            custom_pdf = @(gx, ga, gr, gs)(robust_blurred_uniform_pdf(gx, ga, gr, gs)) .* (1 - minimum_density) + minimum_density;
          case 'constrained_beta'
            custom_pdf = @(gx, ga, gb)(betapdf(gx, ga, gb)) .* (1 - minimum_density) + minimum_density;
          case 'raw'
            % custom_pdf = @(gx, gp)(result_parameters(find(gx < gp, 1, 'last'))) .* (1 - minimum_density) + minimum_density;
            % This is hacky, but it's easier than making an exception to the "result_densities = " line below. Sometime I should make these models classes or at least structs with a uniform interface:
            % custom_pdf = @(gx, varargin)(varargin{find(gx < radial_raw_bin_edges, 1, 'last')}) .* (1 - minimum_density) + minimum_density;
            custom_pdf = @(gx, varargin)subsref(cell2mat(varargin), struct('type', '()', 'subs', {{arrayfun(@(gxx)find(gxx >= radial_raw_bin_edges, 1, 'last'), gx)}})) .* (1 - minimum_density) + minimum_density;
            % bar(radial_raw_bin_edges, result_parameters, 'histc')
            % error('Finish implementing raw!!!')
          otherwise
            error('Unrecognized current_radial_distribution_type "%s"', current_radial_distribution_type)
        end
        result_densities = custom_pdf(given_radial_coordinates(current_valid_region), current_parameters_cell_without_scale{:});
        
      otherwise
        result_densities = pdf(current_radial_distribution_type, given_radial_coordinates(current_valid_region), current_parameters_cell_without_scale{:});
    end
    result_densities(~isfinite(result_densities)) = nan;
  end
  
  
  function [result_image] = fwd_cyto_radial_model_reconstruction(given_parameters, current_radial_distribution_type)
    current_valid_region = approximate_cytoplasm_slice_template_image;
    % current_parameters_cell = num2cell(given_parameters);
    current_parameters_cell_without_scale = num2cell(given_parameters(1:end - 1));
    % current_densities = pdf(current_radial_distribution_type, template_cylindrical_r(current_valid_region), subsref(current_parameters_cell(1:end - 1), struct('type', '{}', 'subs', {{':'}})));
    % Hack values so densities aren't actually computed at 0 or 1:
    template_cylindrical_r = min(template_cylindrical_r, 1 - eps * 1e10);
    template_cylindrical_r = max(template_cylindrical_r, eps * 1e10);
    switch current_radial_distribution_type
      case 'normal_curve'
        curve_fit = cfit(fittype('gauss1'), current_parameters_cell_without_scale{:})
        current_densities = curve_fit(template_cylindrical_r(current_valid_region)).';
        
      case {'truncated_normal', 'blurred_uniform', 'constrained_beta'}
        % Don't duplicate code like this!
        switch current_radial_distribution_type
          case 'truncated_normal'
            % Know want truncated to range [0, 1], so not free params:
            % custom_pdf = @(gx, gm, gs)normpdf(gx, gm, gs) .* ((gx >= 0) & (gx <= 1)) ./ diff(normcdf([0, 1], gm, gs));
            custom_pdf = @(gx, gm, gs)(normpdf(gx, gm, gs) .* ((gx >= 0) & (gx <= 1)) ./ diff(normcdf([0, 1], gm, gs))) .* (1 - minimum_density) + minimum_density;
          case 'blurred_uniform'
            % custom_pdf = @robust_blurred_uniform_pdf;
            custom_pdf = @(gx, ga, gr, gs)(robust_blurred_uniform_pdf(gx, ga, gr, gs)) .* (1 - minimum_density) + minimum_density;
          case 'constrained_beta'
            custom_pdf = @(gx, ga, gb)(betapdf(gx, ga, gb)) .* (1 - minimum_density) + minimum_density;
          otherwise
            error('Unrecognized current_radial_distribution_type "%s"', current_radial_distribution_type)
        end
        current_densities = custom_pdf(template_cylindrical_r(current_valid_region), current_parameters_cell_without_scale{:});
        % error
        
      case 'raw'
        % custom_pdf = @(gx, varargin)(varargin{find(gx < radial_raw_bin_edges, 1, 'last')}) .* (1 - minimum_density) + minimum_density;
        % bar(radial_raw_bin_edges, cell2mat(current_parameters_cell_without_scale), 'histc')
        bar(radial_raw_bin_edges(1:end - 1), cell2mat(current_parameters_cell_without_scale), 'histc')
        for bin_index = 1:length(radial_raw_bin_edges) - 1
          current_densities((template_cylindrical_r(current_valid_region) >= radial_raw_bin_edges(bin_index + 0)) & (template_cylindrical_r(current_valid_region) < radial_raw_bin_edges(bin_index + 1))) = current_parameters_cell_without_scale{bin_index};
        end
        if false
          % Debug info:
          beep
          dbstack
          % pause
          keyboard
        end
        % error('Finish implementing raw!!!')
      otherwise
        current_densities = pdf(current_radial_distribution_type, template_cylindrical_r(current_valid_region), current_parameters_cell_without_scale{:});
        % % Does this normalization even work? Nvm, moved up into normalization parameter:
        % current_cdf_at_furthest = cdf(current_radial_distribution_type, max(template_cylindrical_r(current_valid_region)), current_parameters_cell{1:end - 1});
        % result_image = result_image ./ current_cdf_at_furthest;
    end
    current_densities(~isfinite(current_densities)) = nan;
    result_image = extrapolate_nans(subsasgn(nan(size(template_image)), struct('type', '()', 'subs', {{current_valid_region}}), current_densities)) .* template_image .* given_parameters(end);
    if false && master_script_options.debug_model_reconstruction
      imagesc(reshape_contrast(result_image, -1)), colormap(illustration_colormap_to_use), colorbar, axis equal
      beep
      pause
      % keyboard
    end
  end
  
  
  for model_type_index = 1:number_model_types
    model_type = model_types{model_type_index};
    switch model_type
      case 'standardized_voxels'
        % 2013-08-11 this currently uses 6628 / (71*71*35) = 0.0376 of the space used by all voxels:
        model_representation_functions{model_type_index} = @(x)reshape(x(template_image), 1, []);
        model_reconstruction_functions{model_type_index} = @(x)subsasgn(zeros(size(template_image)), struct('type', '()', 'subs', {{template_image}}), x);
      case 'marginal_x'
        % Take the mean
        % 2013-09-25 how many parameters?:
        model_representation_functions{model_type_index} = @(x)reshape(sum(sum(x(template_occupied_x_slices, :, :) .* template_image(template_occupied_x_slices, :, :), 2), 3) ./ sum(sum(template_image(template_occupied_x_slices, :, :), 2), 3), 1, []);
        model_reconstruction_functions{model_type_index} = @(x)repmat(template_occupied_x_slice_pad_function(x), [0, 1, 1] .* size(template_image) + [1, 0, 0]) .* template_image;
        model_dimensionality_reduction_methods{model_type_index} = 'pca';
        model_dimensionality_reduction_method_parameters{model_type_index} = struct('explained_variance', .95);
      case 'synapse_yz'
        % 2013-09-25 how many parameters?:
        model_representation_functions{model_type_index} = @(x)reshape(x(approximate_synapse_slice_template_image), 1, []);
        % model_reconstruction_functions{model_type_index} = @(x)subsasgn(zeros(size(template_image)), struct('type', '()', 'subs', {{approximate_synapse_slice_template_image}}), x);
        % Extrapolate values for better contrast stretching during visualization:
        model_reconstruction_functions{model_type_index} = @(x)extrapolate_nans(subsasgn(nan(size(template_image)), struct('type', '()', 'subs', {{approximate_synapse_slice_template_image}}), x)) .* template_image;
        model_dimensionality_reduction_methods{model_type_index} = 'pca';
        model_dimensionality_reduction_method_parameters{model_type_index} = struct('explained_variance', .95);
        model_illustration_view{model_type_index} = 'front';
      case 'syn_mar_y'
        % 2013-09-26 how many parameters?:
        % error
        model_representation_functions{model_type_index} = @(x)reshape(subsref(sum(sum(x .* approximate_synapse_slice_template_image, 3), 1) ./ sum(sum(approximate_synapse_slice_template_image, 3), 1), struct('type', '()', 'subs', {{template_occupied_y_slices}})), 1, []);
        % Extrapolate values for better contrast stretching during visualization:
        model_reconstruction_functions{model_type_index} = @(x)extrapolate_nans(subsasgn(nan(size(template_image)), struct('type', '()', 'subs', {{':', template_occupied_y_slices, ':'}}), repmat(x, [1, 0, 1] .* size(template_image) + [0, 1, 0]))) .* template_image;
        model_dimensionality_reduction_methods{model_type_index} = 'pca';
        model_dimensionality_reduction_method_parameters{model_type_index} = struct('explained_variance', .95);
        model_illustration_view{model_type_index} = 'front';
      case 'fwd_cyto'
        % 2013-10-03 how many parameters?:
        model_representation_functions{model_type_index} = @(x)reshape(x(approximate_cytoplasm_slice_template_image), 1, []);
        % Extrapolate values for better contrast stretching during visualization:
        model_reconstruction_functions{model_type_index} = @(x)extrapolate_nans(subsasgn(nan(size(template_image)), struct('type', '()', 'subs', {{approximate_cytoplasm_slice_template_image}}), x)) .* template_image;
        % model_dimensionality_reduction_methods{model_type_index} = 'pca';
        model_dimensionality_reduction_method_parameters{model_type_index} = struct('explained_variance', .95);
        model_illustration_view{model_type_index} = 'front';
      % case 'fwd_cyto_radial'
      case radial_distribution_model_types
        radial_distribution_model_type_index = find(strcmp(model_type, radial_distribution_model_types));
        radial_distribution_type = radial_distribution_types{radial_distribution_model_type_index};
        % Parameters are estimated mean and standard deviation of radial position of intensity, fitted by a distribution of type radial_distribution_type, and the mean intensity:
        model_representation_functions{model_type_index} = @(x)fwd_cyto_radial_model_representation(x, radial_distribution_type);
        % Extrapolate values for better contrast stretching during visualization:
        model_reconstruction_functions{model_type_index} = @(x)fwd_cyto_radial_model_reconstruction(x, radial_distribution_type);
        % model_dimensionality_reduction_methods{model_type_index} = 'pca';
        model_dimensionality_reduction_method_parameters{model_type_index} = struct('explained_variance', .95);
        model_illustration_view{model_type_index} = 'front';
      case 'cylindrical_zernike'
        zernike_options = template_options;
        % zernike_options.highest_order_same_across_slices = true;
        % This produces a high intensity line down the axis for some reason:
        % highest_order = ceil((zernike_options.yr + zernike_options.zr) * .25);
        % Maximum should be 15 parameters per slice:
        % highest_order = 4;
        % Maximum should be 21 parameters per slice:
        % highest_order = 5;
        highest_order = 6;
        % highest_order = 8;
        % Should add X axis downsampling sometime...
        % cylindrical_zernike_x_scale = .25;
        % cylindrical_zernike_x_scale = .5;
        cylindrical_zernike_x_scale = 1;
        zernike_options.xc = (zernike_options.xc - 1) * cylindrical_zernike_x_scale + 1;
        zernike_options.xr = zernike_options.xr * cylindrical_zernike_x_scale;
        zernike_options.imx = ceil(zernike_options.imx * cylindrical_zernike_x_scale);
        % model_should_compute_covariance(model_type_index) = true;
        model_dimensionality_reduction_methods{model_type_index} = 'pca';
        % warning('>>> HACK, model_dimensionality_reduction_methods for %s = none!', model_type), model_dimensionality_reduction_methods{model_type_index} = 'none';
        model_dimensionality_reduction_method_parameters{model_type_index} = struct('explained_variance', .95);
        % model_dimensionality_reduction_method_parameters{model_type_index} = struct('explained_variance', .99);
        if cylindrical_zernike_x_scale == 1
          model_representation_functions{model_type_index} = @(x)reshape(compute_template_pseudo_zernike_moments(x, highest_order, zernike_options), 1, []);
          model_reconstruction_functions{model_type_index} = @(x)reconstruct_template_image_from_pseudo_zernike_moments(x, highest_order, zernike_options);
        else
          model_representation_functions{model_type_index} = @(x)reshape(compute_template_pseudo_zernike_moments(image_resize_nd(x, [cylindrical_zernike_x_scale, 1, 1], 'bilinear'), highest_order, zernike_options), 1, []);
          model_reconstruction_functions{model_type_index} = @(x)image_resize_nd(reconstruct_template_image_from_pseudo_zernike_moments(x, highest_order, zernike_options), [template_options.imx / zernike_options.imx, 1, 1], 'bilinear');
        end
        if false
          slice_highest_orders = get_template_highest_orders(highest_order, options)
          error
        end
        if false
          downsampled_template = image_resize_nd(double(template_image), [.5, 1, 1], 'bilinear');
          downsampled_template_row_marginals = sum(sum(downsampled_template, 2), 3)
          imshow(reshape_contrast(downsampled_template, -1))
          error
        end
      case 'spherical_zernike'
        error('Not yet implemented!')
    end
  end

  if master_script_options.debug_model_reconstruction
    % % Test model_type = 'standardized_voxels' for values inside the template:
    % Test each model_type for values inside the template:
    clf, set(gcf, 'Position', [1, 1, 4, 3] * 200);
    % movegui(gcf, 'center');
    % for model_type_index = 1:number_model_types
    % warning('>>>> HACK, just test fwd_cyto_radial'), for model_type_index = 6
    warning('>>>> HACK, just test fwd_cyto_radial_raw'), for model_type_index = 10
      model_type = model_types{model_type_index}
      for a_index = 1:8
        switch a_index
          case 1
            % Should have zero error:
            a = double(template_image);
          case 2
            % Should have zero error:
            % a = template_image + template_image .* randn(size(template_image));
            % Modified to have non-negative values so can treat as a density:
            a = template_image + template_image .* rand(size(template_image));
          case 3
            % Should have significant error outside the template:
            % a = template_image + randn(size(template_image));
            % Modified to have non-negative values so can treat as a density:
            a = template_image + rand(size(template_image));
          case 4
            % Distance transform:
            a = bwdist(~template_image);
          case 5
            % Distance transform to generate a set of rings based on the distance transform:
            a = bwdist(~template_image);
            % a = contrast_stretch(a);
            a(template_image) = contrast_stretch(a(template_image));
            a = 1 - sawtooth(double(a) * 2 * pi * 2, .5);
            % a = 1 - sawtooth(double(a) * 2 * pi * 4, .5);
            prctile(a(template_image), 0:25:100)
            a(~template_image) = 0;
            a(template_image) = contrast_stretch(a(template_image));
            % prctile(a(template_image), 0:25:100)
          case 6
            % Distance transform to generate a single peripheral ring:
            template_image_extended = template_image;
            template_image_extended(1:template_options.xc, :, :) = template_cylindrical_r(1:template_options.xc, :, :) <= template_options.xr;
            a = bwdist(~template_image_extended) .* template_image;
            % a = contrast_stretch(a);
            a(template_image) = contrast_stretch(a(template_image));
            a = 1 - sawtooth((double(a) + .5) * pi * 2, .5);
            % a = 1 - sawtooth(double(a) * 2 * pi * 4, .5);
            prctile(a(template_image), 0:25:100)
            a(~template_image) = 0;
            a(template_image) = contrast_stretch(a(template_image));
            % prctile(a(template_image), 0:25:100)
          case 7
            % Directly generate a peripheral ring using sawtooth:
            a = template_cylindrical_r .* template_image;
            % a = contrast_stretch(a);
            a(template_image) = contrast_stretch(a(template_image));
            a = 1 - sawtooth((double(a) + .5) * pi * 2, .5);
            % a = 1 - sawtooth(double(a) * 2 * pi * 4, .5);
            prctile(a(template_image), 0:25:100)
            a(~template_image) = 0;
            a(template_image) = contrast_stretch(a(template_image));
            % prctile(a(template_image), 0:25:100)
          case 8
            % Directly generate a peripheral ring using lognormal:
            a = template_cylindrical_r .* template_image;
            a(template_image) = pdf('lognormal', template_cylindrical_r(template_image), log(.5), .5);
            % a = 1 - sawtooth(double(a) * 2 * pi * 4, .5);
            prctile(a(template_image), 0:25:100)
            a(~template_image) = 0;
            a(template_image) = contrast_stretch(a(template_image));
            % prctile(a(template_image), 0:25:100)
        end
        ap = model_representation_functions{model_type_index}(a);
        ar = model_reconstruction_functions{model_type_index}(ap);
        whos a ap ar
        if length(ap) <= 20
          ap
        end
        % ap_percentiles = prctile(ap(:), 0:25:100)
        % ar_percentiles = prctile(ar(:), 0:25:100)
        a_reconstruction_norm = norm(a(template_image) - ar(template_image))
        a_reconstruction_norm_ratio = a_reconstruction_norm / norm(a(template_image))
        image_to_display = [template_crop_function(a); template_crop_function(ar); template_crop_function(a - ar)];
        % imshow(reshape_contrast(image_to_display, -1))
        % imshow(reshape_contrast(image_to_display, 2))
        % imshow(reshape_contrast(image_to_display, -1), 'Border', 'tight')
        % imagesc(reshape_contrast(image_to_display, -1)), axis image
        imagesc(reshape_2d(image_to_display, 2)), axis image
        % colormap(pmkmp(1024, 'CubicL'))
        colormap(illustration_colormap_to_use)
        title({'Original pattern, reconstruction, difference', sprintf('model_type: %s', model_type)}, 'Interpreter', 'none')
        proportion_zero_valued = mean(ap == 0)
        if false && strcmp(model_type, 'cylindrical_zernike')
          ap
        end
        pause
      end
    end
    error('This test finished')
  end

  
  t_cell_info.model_type_info = struct()
  
  t_cell_info.model_type_info.model_types = model_types;
  t_cell_info.model_type_info.number_model_types = number_model_types;
  t_cell_info.model_type_info.radial_distribution_types = radial_distribution_types;
  t_cell_info.model_type_info.radial_distribution_model_types = radial_distribution_model_types;
  t_cell_info.model_type_info.number_radial_distribution_types = number_radial_distribution_types;
  
  t_cell_info.model_type_info.model_representation_functions = model_representation_functions;
  t_cell_info.model_type_info.model_reconstruction_functions = model_reconstruction_functions;
  
  t_cell_info.model_type_info.model_dimensionality_reduction_methods = model_dimensionality_reduction_methods;
  t_cell_info.model_type_info.model_dimensionality_reduction_method_parameters = model_dimensionality_reduction_method_parameters;
  t_cell_info.model_type_info.model_dimensionality_reduction_method_transformations = model_dimensionality_reduction_method_transformations;
  
  t_cell_info.model_type_info.model_illustration_view = model_illustration_view;
  
  t_cell_info.model_type_info.fwd_cyto_radial_model_density = @fwd_cyto_radial_model_density;
  t_cell_info.model_type_info.radial_raw_number_bins = radial_raw_number_bins;
  % t_cell_info.model_type_info.radial_raw_bin_centers = radial_raw_bin_centers;
  t_cell_info.model_type_info.radial_raw_bin_edges = radial_raw_bin_edges;

end
  
  
