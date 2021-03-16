function [ embedding_model ] = embed_distance_matrix( distances_incomplete, options )
%takes an 
%   n x n sparse matrix distances_incomplete
%   options structure 
    if ~exist('options', 'var') | isempty(options)
        options = struct();
    end
    
    options = ml_initparam(options, struct( ...
        'number_images', max(size(distances_incomplete)), ...
        'embedding', struct(), ...
        'compute_triangulation', true, ...
        'positions', [] ...
        ));
    
    embedding_info = [];
    
    %if there isn't already an embedding
    if isempty(options.positions)
        
        options.embedding = ml_initparam(options.embedding, struct( ...
        'method', 'most_complete', ...
        'force_positive_definiteness', true, ... %% this is only for the embed_partial_distance_matrix function
        'explained_variance_threshold', 0.9, ...
        'desired_shape_space_dimensionality', 7, ...
        'truncate_dimensions', 7 ...
        ));
        
        [positions, mass_matrix, embedding_info] = embed_partial_distance_matrix(distances_incomplete, options.embedding);

        if size(positions,2) > options.embedding.truncate_dimensions
            positions = positions(:,1:options.embedding.truncate_dimensions);
        end
    else
        positions = options.positions;
    end
    
    keepinds = ~any(isnan(positions),2);
    
    
    embedding_model = options;
    embedding_model.positions               = positions;
    
    if options.compute_triangulation
        embedding_model.convex_hull             = convhulln( positions(keepinds,:)); %result.convex_hull;
        embedding_model.tessellation            = delaunayn(positions(keepinds,:)); %result.tes;

        keepinds_temp = find(keepinds);

        embedding_model.convex_hull             = keepinds_temp(embedding_model.convex_hull);
        embedding_model.tessellation            = keepinds_temp(embedding_model.tessellation);
    end

    
    
    disp([num2str(sum(keepinds)) filesep num2str(length(keepinds)) ' points embedded.'])
            

%     positions_temp                            = embedding_model.positions;
%     embedding_model.positions               = nan(options.number_images, size(embedding_model.positions,2));
%     embedding_model.positions(keepinds,:)   = positions_temp;
    embedding_model.keepinds                = keepinds;
    
    
    embedding_model.distances_incomplete    = distances_incomplete;

    embedding_model.embedding_info          = embedding_info;
    embedding_model.embedding               = options.embedding;

end

