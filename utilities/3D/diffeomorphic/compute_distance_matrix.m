function [distances, distances_sampled, success] = compute_distance_matrix(...
    image_function, number_shapes, options)
% [result, success] = compute_distance_matrix(image_function, number_shapes, save_filename_prefix, registration_options, convergence_error_function, desired_shape_space_dimensionality, maximum_parallel_cpus)
%
%   Computes the LDDMM distance between each pair of shapes and returns a distance matrix in result.
%   Parallel.
%
%   2012-12-01 tebuck: copied from register_shapes_to_mean_shape.m.
%   2013-04-20 tebuck: added argument  to stop when the number of distances known is the same as the number of degrees of freedom of a distance matrix for points of a specified dimensionality.
%   2013-04-24 tebuck: adding parfor if new argument maximum_parallel_cpus is greater than one.
%   2013-05-13 tebuck: skipped distances are now silent as those messages were mostly noise.
%   2013-08-14 gj:     Removed skipping distance matrix if file already
%                      exists, so we can add images and update the shape
%                      as preprocessed results complete
%   2013-08-18 gj:     Distance matrix fills randomly to less bias reconstruct_distance_matrix
%   2013-08-28 gj:     Added return statement in catch block of
%                      compute_single_distance, and removed parfor loop
%   2013-08-29 gj:     Added "useCurrentResults" to generate a model with
%                      only current results, and fixed bug where .tmp
%                      files were being left behind
% Aug 30, 2013 G. Johnson    Changed they way files are input into the
%                            diffeomorphic model function
% Sept 4, 2013 G. Johnson Improved error handling
% March 23, 2016 xruan Intergrate faster code for distance computing.
% July 13, 2016 xruan add try block for extract images.

options = ml_initparam(options, struct( ...
    'desired_shape_space_dimensionality', inf, ...
    'useCurrentResults',            false, ...
    'matCompletionFunctionString',  'defaultCompletion(distances_sampled, dist_nans)', ...
    'save_location', [], ...
    'distances_computed_before_saving', 1, ...
    'distance_computing_method', 'faster' ...
    ));


desired_shape_space_dimensionality  = options.desired_shape_space_dimensionality;
%     convergence_error_function          = options.convergence_error_function;
useCurrentResults                   = options.useCurrentResults;
matCompletionFunctionString         = options.matCompletionFunctionString;
save_location                       = options.save_location;
distance_computing_method           = options.distance_computing_method;

distances_sparse_file = [save_location filesep 'distances_sparse.mat'];
distances_save_dir = [save_location filesep 'distances'];

if ~exist(distances_save_dir, 'dir')
    mkdir(distances_save_dir)
end


[distances_incomplete, distances_sampled] = loadCalculatedDistances(number_shapes, distances_save_dir, distances_sparse_file, false);

success = true;
isdone = false;

if useCurrentResults
    disp('useCurrentResults flag is true. Using current distance matrix.')
    isdone = true;
end

dist_nans = logical(sparse(number_shapes, number_shapes));

% Either compute distances in series or just load them and check success using this same loop:
dists_computed_tmp = 0;
ind1_temp = nan(options.distances_computed_before_saving, 1);
ind2_temp = nan(options.distances_computed_before_saving, 1);
distances_temp = nan(options.distances_computed_before_saving, 1);

individual_success = true;
while ~isdone
    
% Ulani Qi (uhq) 01/18/18 added pragma below to help eval run when isdeployed
    %#function matrix_completion_random_columns

    [ind1, ind2, isdone] = eval(matCompletionFunctionString);
    
    
    if ~isdone
        filedir = [distances_save_dir filesep 'd' num2str(ind1)];
        
        if ~exist(filedir, 'dir')
            mkdir(filedir)
        end
        
        save_filename = [filedir, filesep, num2str(ind1, '%d'), num2str(ind2, '_%d')];
        
        %         if ~exist([save_filename '.mat'], 'file')
        
        disp(['Computing distance between (' num2str(ind1) ', ' num2str(ind2) ')']);
        [distance, individual_success, ismissing] = compute_single_distance(image_function, options, ind1, ind2, save_filename, distance_computing_method);
        
        %if any images are missing, then we add them to the distance nans
        %so we dont have to spend time attempting to compute that distance
        if any(ismissing)
            if ismissing(1)
                dist_nans(ind1,:) = true;
                dist_nans(:,ind1) = true;
            end
            
            if ismissing(2)
                dist_nans(ind2,:) = true;
                dist_nans(:,ind2) = true;
            end
        end
        
        if isnan(distance)
            dist_nans(ind1, ind2) = true;
            dist_nans(ind2, ind1) = true;
        else
            distances_sampled(ind1, ind2) = 1;
        end
        
        dists_computed_tmp = dists_computed_tmp +1 ;
        ind1_temp(dists_computed_tmp) = ind1;
        ind2_temp(dists_computed_tmp) = ind2;
        if parallel.gpu.GPUDevice.isAvailable && isa(distance,'gpuArray')
            distances_temp(dists_computed_tmp) = gather(distance);
        else
            distances_temp(dists_computed_tmp) = distance;
        end
    end
    
    if isdone || (dists_computed_tmp >= options.distances_computed_before_saving)
        [distances_incomplete, distances_sampled] = addCalculatedDistances(ind1_temp, ind2_temp, distances_temp, distances_sparse_file);
        
        dists_computed_tmp = 0;
        ind1_temp = nan(options.distances_computed_before_saving, 1);
        ind2_temp = nan(options.distances_computed_before_saving, 1);
        distances_temp = nan(options.distances_computed_before_saving, 1);
    end
    % Distance matrix is not complete if one distance is missing:
    success = success && individual_success;
end

%     [distances_incomplete, distances_sampled] = loadCalculatedDistances(number_shapes, distances_save_dir, distances_sparse_file, true);

distances_incomplete = full(distances_incomplete) + full(distances_incomplete)';
distances_incomplete((full(distances_sampled) + full(distances_sampled)')==0) = NaN;

distances = distances_incomplete;
end

function [distance, success, ismissing] = compute_single_distance(image_function, registration_options, shape_index, shape_index2, save_filename, distance_computing_method)

convergence_error_function = @(s, t)0;

% Only compute the first (desired_shape_space_dimensionality + 2) columns:
%   if shape_index > (desired_shape_space_dimensionality + 2)
%     % continue
%     distance = nan;
%     % success = false;
%     success = true;
%     return
%   end

% Check if work has been done
imind2 = shape_index;
imind1 = shape_index2;

[can_start, final_name, final_exists, tempname] = chunk_start(save_filename);

ismissing = false(1,2);

if ~can_start
    %     stack = dbstack();
    % warning('Skipping "%s" at line %d in file %s.', save_filename, stack(1).line, stack(1).file)
    % fprintf('compute_distance_matrix: Skipping "%s" at line %d in file %s.\n', save_filename, stack(1).line, stack(1).file)
    if final_exists
        success = true;
        r = [];
        
        try
            load([save_filename '.mat'], 'r')
            %       % distances(shape_index, shape_index2) = r.total_distance;
            distance = r.total_distance;
            
        catch
            distance = NaN;
            success = false;
        end
        return;
    else
        success = false;
        
        warning([save_filename ' being computed elsewhere.'])
        fprintf('compute_distance_matrix: Being computed elsewhere.\n')
        
        % distances(shape_index, shape_index2) = nan;
        distance = NaN;
        return;
    end
else
    
    fail = false;
    % add try block for failure case.
    try
        source = [];
        target = [];
        [source] = image_function(shape_index);
        [target] = image_function(shape_index2);
    catch
        fail = true;
    end
    
    ismissing = [isempty(source), isempty(target)];

    if fail || any(ismissing)
        shapes = [shape_index shape_index2];
        
        if isempty(source)
            disp(['Could not retrieve source image(s) ' num2str(shapes(ismissing)) '.'])
            pause(1)
        else
            disp(['Could not retrieve target image(s) ' num2str(shapes(ismissing)) '.'])
            pause(1)
        end
        distance = NaN;
        success = false;
    else
        fprintf('compute_distance_matrix: Computing %s\n', save_filename)
        
        % Perform work
        tic
        registration_options.convergence_registration_error = ...
            convergence_error_function(source, target);
        registration_options.single_sided = false;
        registration_options.just_compute_distance = true;
        % registration_options.verbose = 0;
        % registration_options
        if strcmp(distance_computing_method, 'faster')
            r = Greedy3D_lambda_pre(...
                source, target, 1, registration_options)
            
        else
            r = Greedy3D_lambda_pre_compressed(...
                source, target, 1, registration_options)
        end
        % Save work
        r.time_elapsed = toc;
        
        save([save_filename '.mat'], 'r')
        
        distance = r.total_distance;
        success = true;
    end
    
    chunk_finish(tempname)
end
end

function [ind1, ind2, isdone] = defaultCompletion(distances_sparse_file, dist_nans)

inds = find(triu(~(distances_sparse_file| dist_nans)));

if isempty(inds)
    isdone = true;
    ind1 = [];
    ind2 = [];
    return
end

ind = inds(1);
isdone = false;
[ind1, ind2] = ind2sub(size(distances_sparse_file), ind);
end

function [distances_incomplete, distances_sampled] = addCalculatedDistances(ind1, ind2, distance, distances_sparse_file)
canstart = false;

while ~canstart
    disp('Waiting for distance file lock')
    [canstart, ~, ~, temp_name] = chunk_start([distances_sparse_file '_lock'], '');
    if ~canstart
        pause(1);
    end
end

load(distances_sparse_file)

for i = 1:length(distance)
    if ~isnan(distance)
        distances_incomplete(ind1(i), ind2(i)) = distance(i);
        %     distances_incomplete(ind2, ind1) = distance;
        
        distances_sampled(ind1(i), ind2(i)) = 1;
    end
end
%     distances_sampled(ind2, ind1) = 1;

disp('Writing updated distance file');
save([distances_sparse_file '_bak'], 'distances_incomplete', 'distances_sampled', '-v7.3')
copyfile([distances_sparse_file '_bak'], distances_sparse_file);

disp('Releasing distance file')
chunk_finish(temp_name);
end

function [distances_incomplete, distances_sampled] = loadCalculatedDistances(number_shapes, distances_dir, distances_sparse_file, check_missing)
%if there is not a distances file and there is a lock
%wait for it to show up (because someone else is assembling it)
%end

%if there is a distances file
%load it
%elseif there is a not a distances file
%lock it
%assemble it
%end

canstart = false;

while ~canstart
    disp('Waiting for distance file lock')
    [canstart, ~, ~, temp_name] = chunk_start([distances_sparse_file '_lock'], '');
    if ~canstart
        pause(1);
    end
end

file_exists = exist(distances_sparse_file, 'file');

if file_exists
    disp('Loading existing distances file')
    load(distances_sparse_file)
    
    distances_sampled = (distances_sampled + distances_sampled') ~= 0;
else
    disp('Constructing distance file from current results')
    
    distances_incomplete = sparse(number_shapes, number_shapes);
    
    distances_sampled = sparse(number_shapes, number_shapes);
    distances_sampled(logical(eye(size(distances_incomplete)))) = 1;
end

if ~file_exists || check_missing
    
    files = dir([distances_dir '*' filesep '*.mat']);
    files = {files.name};
    tokens = regexp(files, '[^0-9]', 'split');
    tokens = vertcat(tokens{:});
    
    if ~isempty(tokens)
        
        xy = str2double(tokens(:,1:2));
        x = xy(:,1);
        y = xy(:,2);
        
        
        for i = 1:size(xy(:,1))
            disp(['Loading ' num2str(i) filesep num2str(size(xy(:,1)))])
            if distances_sampled(x(i), y(i)) == 0
                save_filename = [distances_dir filesep 'd' num2str(ind1) filesep num2str(x(i)) '_' num2str(y(i)) '.mat'];
                
                try
                    load(save_filename, 'r')
                    
                    distances_incomplete(x(i), y(i)) = r.total_distance;
                    distances_sampled(x(i), y(i)) = 1;
                catch
                    warning(['Could not read ' save_filename]);
                end
            end
        end
    end
    save(distances_sparse_file, 'distances_incomplete', 'distances_sampled')
    
end

chunk_finish(temp_name);
end
