function [result_structure] = reduce_model_dimensionality(options)
  % Given a function batch_function that gives a subset of available data, the number of subsets number_batches, and a dimensionality reduction method method, return reduction and reconstruction functions for that data and method. Assumes that batches are from different models for at least the 'pca' method so that at least some models can later be compared with statistical tests that have sample size restrictions, e.g., MANOVA (need more data than features, IIRC).
  %
  % Tests:
  %
  % s = RandStream('mt19937ar', 'Seed', 597075);
  % RandStream.setDefaultStream(s);
  % data_batches = arrayfun(@(given_index)randn(randi([50, 100]), 5) * randn(5, 10), (1:5)', 'UniformOutput', false)
  % options = struct();
  % options.number_batches = length(data_batches);
  % options.batch_function = @(given_index)data_batches{given_index};
  % reduced_model_functions = reduce_model_dimensionality(options)
  % all_data = cell2mat(data_batches);
  % mean_mean_absolute_value = mean(abs(all_data(:)))
  % mean_mean_reconstruction_absolute_error = mean(mean(abs(reduced_model_functions.reconstruction_function(reduced_model_functions.representation_function(all_data)) - all_data)))
  % mean_standard_deviation = mean(std(all_data))
  % reconstruction_mean_standard_deviation = mean(std(reduced_model_functions.reconstruction_function(reduced_model_functions.representation_function(all_data)) - all_data))
  %
  % 2013-08-16 tebuck: Copied from master_script.m.
  % % Dependencies:
  % % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
  % % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh

  
  default_options = struct();
  default_options.number_batches = 0;
  default_options.batch_function = @(given_index)[];
  default_options.method = 'pca';
  % Proportion of models to be able to test if too many principal components are needed to explain options.method_parameters.pca.explained_variance of the data to successfully test all models:
  default_options.testable_proportion = .95;
  warning('>>>> HACK, default_options.testable_proportion = .05'), default_options.testable_proportion = .05;
  default_options.save_filename = [];
  default_options.verbosity = 0;

  default_options.method_parameters = struct();
  default_options.method_parameters.none = struct();
  default_options.method_parameters.pca = struct();
  default_options.method_parameters.pca.explained_variance = .95;
  % default_options.method_parameters.pca.explained_variance = .99;
  
  if ~exist('options', 'var')
    options = default_options; 
  else
    options = process_options_structure(default_options, options);
    % options = process_options_structure(default_options, options, [], true);
  end
  % options

  method_need_full_data.none = false;
  % method_need_full_data.pca = false;
  warning('>>>> HACK, method_need_full_data.pca = true!'), method_need_full_data.pca = true;
  
  methods = {'none', 'pca'};
  % methods = fieldnames(method_need_full_data);
  if ~ismember(options.method, methods)
    methods_with_commas = cell2mat(reshape([methods; repmat({', '}, 1, length(methods) - 1), {''}], 1, []));
    error(['Unknown method "%s", choose from {%s}'], options.method, methods_with_commas)
  end

  
  % Compute dimensionality reduction transforms:
  
  transform_parameters = [];
  
  representation_function = [];
  reconstruction_function = [];
  
  % options
  % strcmp(options.method, 'none')
  if strcmp(options.method, 'none')
    representation_function = @(given_parameters)given_parameters;
    reconstruction_function = @(given_parameters)given_parameters;
    % warning('options.method = %s!', options.method)
  else

    % Read all files sequentially:
    need_full_data = method_need_full_data.(options.method);
    if need_full_data
      current_data = [];
    else
      current_number_data = 0;
      switch options.method
        case 'pca'
          current_mean = [];
          current_covariance = [];
      end
    end
    current_minimum_samples = inf;
    current_sample_sizes = [];
    current_number_features = [];
    for batch_index = 1:options.number_batches
      % batch_index
      current_batch = options.batch_function(batch_index);
      % whos current_batch
      if need_full_data
        if isempty(current_data)
          current_data = current_batch;
        else
          current_data = [current_data; current_batch];
        end
      else
        if batch_index == 1
          switch options.method
            case 'pca'
              current_mean = zeros(1, size(current_batch, 2));
              % whos current_mean
              current_covariance = zeros(length(current_mean) .* [1, 1]);
          end
        end
        for datum_index = 1:size(current_batch, 1)
          current_parameters = current_batch(datum_index, :);
          current_number_data = current_number_data + 1;
          switch options.method
            case 'pca'
              % Updated mean:
              preceding_mean = current_mean;
              % whos current_mean
              current_mean = update_running_mean(preceding_mean, current_parameters, current_number_data);
              % whos current_mean
              % Updated covariance:
              preceding_covariance = current_covariance;
              current_covariance = update_running_covariance(preceding_mean, current_mean, preceding_covariance, current_parameters, current_number_data);
            end
        end
        if batch_index == options.number_batches
          switch options.method
            case 'pca'
              current_covariance = finalize_running_covariance(current_covariance, current_number_data);
          end
        end
        % current_number_data
        if options.verbosity > 0
          fprintf('batch_index = %d, current_number_data = %d%s\n', batch_index, current_number_data, repmat('.', 1, 3 * (batch_index < options.number_batches)))
        end
      end
      
      current_sample_sizes(end + 1) = size(current_batch, 1);
      current_number_features = size(current_batch, 2);
    end
    % current_sample_sizes
    current_minimum_samples = min(current_sample_sizes);
    % current_minimum_samples
    if length(current_data) >= 2
      switch options.method
        case 'pca'
          % Proportion of models to be able to test if too many principal components are needed to explain options.method_parameters.pca.explained_variance of the data to successfully test all models:
          if options.testable_proportion > 0
            testable_proportion_number_samples = floor(quantile(current_sample_sizes, 1 - options.testable_proportion));
          else
            % testable_proportion_number_samples = max(current_sample_sizes);
            % Don't reduce dimensionality based on number of samples:
            testable_proportion_number_samples = current_number_features;
          end
          
          if need_full_data
            current_mean = mean(current_data);
            current_covariance = cov(current_data);
            % [~, current_pca_info] = pca(current_data, options.method_parameters.pca.explained_variance);
            % if current_pca_info.num_pcs >= testable_proportion_number_samples
              % desired_explained_variance = sum(current_pca_info.latent(1:floor(testable_proportion_number_samples))) ./ sum(current_pca_info.latent)
              % [~, current_pca_info] = pca(current_data, desired_explained_variance);
            % end
          else
          end
          % Use pcacov on the correlation matrix, which can be computed incrementally, in case the data gets too big at some point or we want to run multiple jobs in limited memory:
          % Reference: <http://www.mathworks.com/help/stats/pcacov.html>
          current_standard_deviations = sqrt(diag(current_covariance));
          current_correlation = current_covariance ./ (current_standard_deviations * current_standard_deviations');
          current_pca_info = struct();
          [current_pca_info.coeff, current_pca_info.latent] = pcacov(current_correlation);
          current_pca_info.num_pcs = min(find(cumsum(current_pca_info.latent) ./ sum(current_pca_info.latent) >= options.method_parameters.pca.explained_variance));
          if current_pca_info.num_pcs >= testable_proportion_number_samples
            current_pca_info.num_pcs = floor(testable_proportion_number_samples);
          end
          current_pca_info.pcs_prop_var = sum(current_pca_info.latent(1:current_pca_info.num_pcs)) ./ sum(current_pca_info.latent);
          % current_pca_info
          
          % Modify representation_function and reconstruction_function:
          % whos
          % options
          % representation_function = @(x)(x - current_mean) * current_pca_info.coeff(:, 1:current_pca_info.num_pcs);
          % reconstruction_function = @(x)x * current_pca_info.coeff(:, 1:current_pca_info.num_pcs).' + current_mean;
          sparse_diagonal_standardization_matrix = spdiags(reshape(1 ./ current_standard_deviations, [], 1), 0, length(current_standard_deviations), length(current_standard_deviations));
          sparse_diagonal_unstandardization_matrix = spdiags(reshape(current_standard_deviations, [], 1), 0, length(current_standard_deviations), length(current_standard_deviations));
          % whos current_mean
          representation_function = @(x)(x - repmat(current_mean, size(x, 1), 1)) * sparse_diagonal_standardization_matrix * current_pca_info.coeff(:, 1:current_pca_info.num_pcs);
          reconstruction_function = @(x)(x * current_pca_info.coeff(:, 1:current_pca_info.num_pcs).') * sparse_diagonal_unstandardization_matrix + repmat(current_mean, size(x, 1), 1);
          
          transform_parameters = current_pca_info;
          
          % error
          
        otherwise
          error('Not yet implemented')
      end
    end

  end
  
  
  result_structure = struct();
  result_structure.representation_function = representation_function;
  result_structure.reconstruction_function = reconstruction_function;
  % result_structure.pca_info = current_pca_info;
  result_structure.transform_parameters = transform_parameters;
  
end
  
  
