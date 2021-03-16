function [r_structure] = Greedy3D_integration(source, target, options)
% Internal function for Greedy3D_lambda_pre_compressed.
%
% This function integrates the differential equations for LDDMM to assist Greedy3D_lambda_pre_compressed.m.
%
% Arguments:
% compressed_source, instance of WindowedImage (if the source and target use compression, this method will do the same to save memory), the moving image (registration) or first/left/source image (interpolation)
% compressed_target, instance of WindowedImage, the fixed image (registration) or second/right/target image (interpolation)
% options, structure, see Greedy3D_compressed_default_options.m for descriptions of options
%
%
% 2011-06-05 15:40 tebuck: Created.
% 2011-06-08 13:23 tebuck: Adding adaptive step size:
% Fehlberg 1970 "Some experimental results concerning the error
% propagation in Runge-Kutta type integration formulas"
% 2011-06-13 14:53 tebuck: Adding outputs of all times and
% distances at those times. Need to add a convergence criterion.
% 2011-06-14 02:22 tebuck: Change in total distance added per
% unit time being less than convergence_tolerance is considered
% convergence.
% 2011-06-16 21:02 tebuck: Due to the relative error tolerance,
% the first step is always accepted. Otherwise, with the distance
% starting at zero, it's amazing that integration ever
% proceeded...
% 2011-06-17 13:59 tebuck: Removing relative error tolerance
% altogether. It's useless in this application, and besides I
% have it saved in revision 208. Step size changes from
% Shacham and Brauner 2007, "Preventing oscillatory behavior in
% error control for ODEs."
% 2011-07-04 18:11 tebuck: Trying to prevent some interpolation
% artifacts by clamping deformations to remain within the
% image. Currently, the distance traveled by the deformation is
% not be computed, so distance is now an upper bound. I plan to
% fix this if the artifacts go away.
% 2011-07-15 14:05 tebuck: New option
% maximum_registration_error_failures quits after that many
% failures. Zero means ignore such failures.
% 2011-07-21 10:23 tebuck: added
% registration_error_failure_step_scale to try to recover when
% registration error starts to inrease.
% 2012-09-20 tebuck: modifying to use WindowedImage for optional in-memory compression.
% 2012-10-24 yueyu: moved the Greedy3D_integration_compressed folder into temp
% 2013-02-09 tebuck: adding known_distance and stop_source_distance options to stop early during interpolation of the distance is already (approximately) known.
% 2013-02-17 tebuck: Changing the default options to those that are currently most useful or are being used in CellOrganizer.
% 2013-02-18 tebuck: Greedy3D_integration_compressed is only called by Greedy3D_lambda_pre_compressed, so simplifying arguments to Greedy3D_integration_compressed.

% prior_warning_state = warning('query', 'MATLAB:DELETE:FileNotFound');
prior_warning_state = warning('off', 'MATLAB:DELETE:FileNotFound');

start_wall_time = tic;
start_cpu_time = cputime;

% Prevent collision by jobs starting simultaneously (this should be an option that lets this call a function to get the job ID so other queue systems can be used):
n = getenv('PBS_JOBID');
n = regexprep(n, '[\r\n]', '');
base_filename = [options.tempparent, filesep, mfilename, filesep, datestr(now(), 'yyyymmddHHMMSSFFF'), '_', n];
warning('off', 'MATLAB:MKDIR:DirectoryExists');
mkdir(fileparts(base_filename))
% warning('on', 'MATLAB:MKDIR:DirectoryExists');


% How large the subimages to process will be:
window_size = options.window_size;
source_size = size(source);
if ndims(source) == 2
    source_size_xyz = [source_size([2, 1]), 1];
    M = source_size(1);
    N = source_size(2);
    P = 1;
elseif ndims(source) == 3
    source_size_xyz = source_size([2, 1, 3]);
    M = source_size(1);
    N = source_size(2);
    P = source_size(3);
end

    function value = compute_registration_error(image1, image2)
        value = sum((image1(:) - image2(:)).^2);
    end

% Default option values and options processing (these are now set in Greedy3D_compressed_default_options.m):

% % Runge-Kutta integrator parameters (Bogacki-Shampine method):
% % Step size fractions:
% default_options.integrator_c = [0, 1/2, 3/4, 1];
% % Stage coefficients for the two final estimates (b is the integration formula, b2 the error/comparison formula):
% default_options.integrator_b = [2/9, 1/3, 4/9, 0];
% default_options.integrator_b2 = [7/24, 1/4, 1/3, 1/8];
% % Per-stage coefficients for the linear combination of prior stages:
% default_options.integrator_a = zeros(4, 3);
% default_options.integrator_a(2, 1) = 1/2;
% default_options.integrator_a(3, 1:2) = [0, 3/4];
% default_options.integrator_a(4, 1:3) = [2/9, 1/3, 4/9];
% % Tells us whether the last stage of one iteration is the same as
% % the first stage of the next iteration:
% default_options.integrator_first_same_as_last = true;
% % Order of error in integration results (is this right?):
% default_options.integrator_error_order = 3;

% % Initial integration step size:
% default_options.step_size = 0.01;

% % Convergence criteria:
% % Threshold of distance moved per unit time below which to declare convergence:
% default_options.convergence_tolerance = 0;
% % Threshold of sum of squared errors between the deformed source and target pixels:
% default_options.convergence_registration_error = 0;
% % Threshold of sum of absolute errors between the deformed source and target pixels:
% default_options.convergence_absolute_error = 0;

% % known_distance says how far the shapes are known/assumed to be (negative values disable this function):
% default_options.known_distance = -1;
% % stop_source_distance is how far the source can move before stopping, and (known_distance - stop_source_distance) is how far the target can move before stopping (negative values disable this function):
% default_options.stop_source_distance = -1;

% % Step size control criteria:
% % Maximum LDDMM distance error (between the two integrations, if applicable) change allowed per iteration:
% default_options.absolute_distance_tolerance = 0;
% % Maximum deformation error in any coordinate for any pixel allowed per iteration:
% default_options.absolute_deformation_tolerance = 0;

% % Switches for saving intermediate results:
% default_options.save_intermediates = false;
% default_options.save_intermediate_images = false;
% default_options.save_stages = true;

% % No interpolation need be returned if just computing distance between images:
% default_options.just_compute_distance = false;
% % True if registration, false if interpolation:
% default_options.single_sided = false;
% % Convergence should not be dependent on number of iterations, so set this to a large number:
% default_options.maximum_iterations = 10000;
% % If the image space is not toroidal/periodic, don't use periodic convolution so, e.g., the left and right sides of images do not interact (if this is false):
% default_options.periodic_space = false;
% % Old option that quits if the registration error goes up too many times (disabled if 0):
% default_options.maximum_registration_error_failures = 0;
% default_options.registration_error_failure_step_scale = 0;
% default_options.maximum_deformation_per_step = ones(1, 3) * .5;
% default_options.verbose = 1;
% % 1x2 double array, the time over which to integrate the differential equations that deform the input images if convergence conditions are not met (set to [0, inf] by default because now there are convergence criteria; by the way, max_time has been ignored for a long time):
% default_options.time_span = [0, inf];
default_options = Greedy3D_compressed_default_options();
if ~exist('options', 'var')
    options = default_options;
else
    options = process_options_structure(default_options, options);
end

current_step_size = options.step_size;
previous_step_size = current_step_size;

% Make maximum_deformation_per_step apply to a 3D image:
if numel(options.maximum_deformation_per_step) == 1
    options.maximum_deformation_per_step = options.maximum_deformation_per_step .* [1, 1, 1];
elseif numel(options.maximum_deformation_per_step) ~= 3
    error('numel(options.maximum_deformation_per_step) ~= 3')
end

if options.verbose > 1
    options.integrator_a, options.integrator_b, options.integrator_b2, options.integrator_c, options.integrator_first_same_as_last, options.integrator_error_order
end


window_mode = 'circular';
if (~options.periodic_space)
    window_mode = 'replicate';
end

current_source_distance = 0;
current_target_distance = 0;

current_source = source;
[X_0,Y_0,Z_0] = meshgrid(1:N, 1:M, 1:P);
% Keep image and window size and whether compressed the same as compressed_source:
%   X_0 = compressed_source.set_data(X_0);
%   Y_0 = compressed_source.set_data(Y_0);
%   Z_0 = compressed_source.set_data(Z_0);
identity_map = {X_0, Y_0, Z_0};
current_source_distance = 0;
current_target_distance = 0;
current_source_deformation = identity_map;
current_target_deformation = [];
current_target = [];
current_target = target;
if ~options.single_sided
    current_target_deformation = current_source_deformation;
end

converged = false;
last_reject_iteration = -Inf;
last_rejected_step_size = current_step_size;
number_rejected_steps = 0;

current_time = options.time_span(1);
time_points = current_time;
source_distance_array = 0;
target_distance_array = 0;

registration_error_array = compute_registration_error(current_source, current_target);
absolute_error_array = sum(abs(current_source(:)-current_target(:)));

iteration_index = 1;
iteration_filenames = cell(1, 0);

if options.just_compute_distance || ~options.save_stages
    iteration_stages = {};
end
if ~options.just_compute_distance && ~options.save_intermediates
    iteration_steps = {};
end

if options.just_compute_distance
elseif options.save_intermediates
    current_iteration_filename = [...
        base_filename ...
        num2str(0, '_iter%05d')...
        ];
    
    total_wall_time = toc(start_wall_time);
    total_cpu_time = cputime - start_cpu_time;
    if options.just_compute_distance
    else
        save([current_iteration_filename '.mat'], ...
            'current_source', 'current_target', ...
            'current_source_deformation', 'current_target_deformation', ...
            'current_source_distance', 'current_target_distance', 'number_rejected_steps', ...
            'total_cpu_time', 'total_wall_time' );
    end
    iteration_filenames = [iteration_filenames, ...
        {current_iteration_filename}];
else
    temp = struct();
    temp.current_source = current_source;
    temp.current_target = current_target;
    temp.current_source_deformation = current_source_deformation;
    temp.current_target_deformation = current_target_deformation;
    temp.current_source_distance = current_source_distance;
    temp.current_target_distance = current_target_distance;
    iteration_steps = [iteration_steps, temp];
end

% Print information headers:
if options.verbose > 0
    integer_format = '% 8d';
    float_format = '% 1.1e';
    info_headers = {'Iter'};
    info_formats = {integer_format};
    if options.convergence_registration_error > 0
        info_headers{end + 1} = 'RegErr';
        info_formats{end + 1} = float_format;
    end
    if options.convergence_absolute_error > 0 || options.convergence_absolute_error_difference > 0
        info_headers{end + 1} = 'AbsErr';
        info_formats{end + 1} = float_format;
    end
    if options.convergence_tolerance > 0
        info_headers{end + 1} = 'DisSpd';
        info_formats{end + 1} = float_format;
    end
    info_headers = [info_headers, {'Time', 'Step'}];
    info_formats = [info_formats, {float_format, float_format}];
    for info_index = 1:length(info_headers)
        fprintf([info_headers{info_index}, repmat(' ', 1, 9 - length(info_headers{info_index}))]);
    end
    fprintf('\n');
    diary_state = get(0, 'Diary'); diary('off'); diary(diary_state);
end

while ~converged && current_time < options.time_span(2)
    
    
    % Print information:
    if options.verbose > 0
        for info_index = 1:length(info_headers)
            info_value = nan;
            switch info_headers{info_index}
                case 'Iter'
                    info_value = iteration_index;
                case 'RegErr'
                    info_value = registration_error_array(end);
                case 'AbsErr'
                    info_value = absolute_error_array(end);
                case 'Time'
                    info_value = current_time;
                case 'Step'
                    info_value = current_step_size;
                case 'DisSpd'
                    if iteration_index == 1
                        info_value = nan;
                    else
                        source_distance_difference = source_distance_array(end) - source_distance_array(end - 1);
                        target_distance_difference = 0;
                        if (~options.single_sided)
                            target_distance_difference = target_distance_array(end) - target_distance_array(end - 1);
                        end
                        info_value = (source_distance_difference + target_distance_difference) ./ (time_points(end) - time_points(end-1)) ./ (source_distance_array(end) + target_distance_array(end));
                    end
            end
            fprintf([info_formats{info_index}, ' '], info_value);
        end
        fprintf('\n');
        diary_state = get(0, 'Diary'); diary('off'); diary(diary_state);
    end
    
    
    % Here we will go through each row of options.integrator_a and compute the
    % stages. We will write intermediate results to disk to save
    % memory.
    stage_filenames = cell(size(options.integrator_a, 1), 1);
    number_stages = size(options.integrator_a, 1);
    % Delete unnecessary stages and steps:
    if iteration_index > 2
        if options.just_compute_distance || ~options.save_stages
            iteration_stages{iteration_index - 2} = [];
        else
            for prior_stage_index = 1:number_stages
                prior_stage_filename = [...
                    base_filename ...
                    num2str(iteration_index - 2, '_iter%05d')...
                    num2str(prior_stage_index, '_stage%05d')...
                    ];
                delete([prior_stage_filename '.mat'])
            end
        end
        if ~options.just_compute_distance && ~options.save_intermediates
            iteration_steps{iteration_index - 2} = [];
            % The first entry is the identity map added before the integration loop:
            iteration_steps{iteration_index - 1} = [];
        end
    end
    for stage_index = 1:number_stages
        start_time = current_time + options.integrator_c(stage_index) * current_step_size;
        start_source_deformation = current_source_deformation;
        start_target_deformation = current_target_deformation;
        start_source_distance = current_source_distance;
        start_target_distance = current_target_distance;
        
        start_source = current_source;
        start_target = current_target;
        
        % Stages are derivatives at the iteration's start point plus
        % linear combinations of previous stages:
        for prior_stage_index = 1:(stage_index - 1)
            prior_stage_filename = [...
                base_filename ...
                num2str(iteration_index, '_iter%05d')...
                num2str(prior_stage_index, '_stage%05d')...
                ];
            if options.just_compute_distance || ~options.save_stages
                temp = iteration_stages{iteration_index}(prior_stage_index);
            else
                temp = load([prior_stage_filename '.mat']);
            end
            prior_source_velocity = temp.d_source;
            prior_target_velocity = temp.d_target;
            prior_source_distance = temp.d_source_distance;
            prior_target_distance = temp.d_target_distance;
            
            % Deformation deformation, not deformation addition:
            
            % The derivatives are actually deformations, so rather than
            % scaling them from zero we add a scaled version to an
            % identity map:
            for dim_ind = 1:3
                prior_source_velocity{dim_ind} = identity_map{dim_ind} + prior_source_velocity{dim_ind} .* options.integrator_a(stage_index, prior_stage_index);
                if (~options.single_sided)
                    prior_target_velocity{dim_ind} = identity_map{dim_ind} + prior_target_velocity{dim_ind} .* options.integrator_a(stage_index, prior_stage_index);
                end
                
            end
            % "Add" the deformation derivative to this stage's start
            % point by interpolation:
            for dim_ind = 1:3
                start_source_deformation{dim_ind} = image_interp3(...
                    start_source_deformation{dim_ind}, prior_source_velocity ...
                    , '*linear', 0, -1, window_mode);
                if (~options.single_sided)
                    start_target_deformation{dim_ind} = image_interp3(...
                        start_target_deformation{dim_ind}, prior_target_velocity ...
                        , '*linear', 0, -1, window_mode);
                end
                
                % Prevent deformations from reaching outside the image:
                start_source_deformation{dim_ind} = max(...
                    start_source_deformation{dim_ind}, 1);
                start_source_deformation{dim_ind} = min(...
                    start_source_deformation{dim_ind}, source_size_xyz(dim_ind));
                if (~options.single_sided)
                    start_target_deformation{dim_ind} = max(...
                        start_target_deformation{dim_ind}, 1);
                    start_target_deformation{dim_ind} = min(...
                        start_target_deformation{dim_ind}, source_size_xyz(dim_ind));
                end
            end
            start_source_distance = start_source_distance + ...
                prior_source_distance * ...
                options.integrator_a(stage_index, prior_stage_index);
            if (~options.single_sided)
                start_target_distance = start_target_distance + ...
                    prior_target_distance * ...
                    options.integrator_a(stage_index, prior_stage_index);
            end
            clear temp prior_source_velocity prior_target_velocity
        end
        
        
        start_source = image_interp3(...
            source, start_source_deformation ...
            , '*linear', 0, -1, window_mode);
        if (~options.single_sided)
            start_target = image_interp3(...
                target, start_target_deformation...
                , '*linear', 0, -1, window_mode);
        end
        d_source = [];
        d_target = [];
        d_source_distance = [];
        d_target_distance = [];
        
        % "First same as last" means that the last stage of the
        % previous iteration is the first stage of this iteration
        % (though due to my implementation I have to rescale the
        % stage to the proper step size below):
        if (options.integrator_first_same_as_last && stage_index == 1 && ...
                iteration_index > 1)
            last_stage_filename = [...
                base_filename ...
                num2str(iteration_index - 1, '_iter%05d')...
                num2str(number_stages, '_stage%05d')...
                ];
            
            if options.just_compute_distance || ~options.save_stages
                d_source = iteration_stages{iteration_index - 1}(number_stages).d_source;
                d_target = iteration_stages{iteration_index - 1}(number_stages).d_target;
                d_source_distance = iteration_stages{iteration_index - 1}(number_stages).d_source_distance;
                d_target_distance = iteration_stages{iteration_index - 1}(number_stages).d_target_distance;
            else
                load([last_stage_filename '.mat'])
            end
            
            for dim_ind = 1:3
                d_source{dim_ind} = d_source{dim_ind} .* (current_step_size / previous_step_size);
                if (~options.single_sided)
                    d_target{dim_ind} = d_target{dim_ind} .* (current_step_size / previous_step_size);
                end
            end
            d_source_distance = d_source_distance * current_step_size / (previous_step_size);
            d_target_distance = d_target_distance * current_step_size / (previous_step_size);
        else
            % Not first same as last, so compute the new derivative:
            [d_source, d_target, d_source_distance, d_target_distance] = ...
                Greedy3D_single_step(start_source, start_target, options);
            for dim_ind = 1:3
                d_source{dim_ind} = d_source{dim_ind} .* current_step_size;
                if (~options.single_sided)
                    d_target{dim_ind} = d_target{dim_ind} .* current_step_size;
                end
            end
            d_source_distance = d_source_distance * current_step_size;
            d_target_distance = d_target_distance * current_step_size;
        end
        
        % Save this stage:
        current_stage_filename = [...
            base_filename ...
            num2str(iteration_index, '_iter%05d')...
            num2str(stage_index, '_stage%05d')...
            ];
        if options.just_compute_distance || ~options.save_stages
            temp = struct();
            temp.d_source = d_source;
            temp.d_target = d_target;
            temp.d_source_distance = d_source_distance;
            temp.d_target_distance = d_target_distance;
            iteration_stages{iteration_index}(stage_index) = temp;
        else
            save([current_stage_filename '.mat'], ...
                'd_source', 'd_target', ...
                'd_source_distance', 'd_target_distance')
        end
        
        clear d_source d_target d_source_distance d_target_distance
    end
    
    % Compute the final new points:
    b_source_deformation = current_source_deformation;
    b_target_deformation = current_target_deformation;
    b2_source_deformation = current_source_deformation;
    b2_target_deformation = current_target_deformation;
    b_source_distance = current_source_distance;
    b_target_distance = current_target_distance;
    b2_source_distance = current_source_distance;
    b2_target_distance = current_target_distance;
    for stage_index = 1:number_stages
        current_stage_filename = [...
            base_filename ...
            num2str(iteration_index, '_iter%05d')...
            num2str(stage_index, '_stage%05d')...
            ];
        
        if options.just_compute_distance || ~options.save_stages
            d_source = iteration_stages{iteration_index}(stage_index).d_source;
            d_target = iteration_stages{iteration_index}(stage_index).d_target;
            d_source_distance = iteration_stages{iteration_index}(stage_index).d_source_distance;
            d_target_distance = iteration_stages{iteration_index}(stage_index).d_target_distance;
        else
            load([current_stage_filename '.mat'])
        end
        % Deformation deformation:
        d_source2 = d_source;
        d_target2 = d_target;
        for dim_ind = 1:3
            d_source{dim_ind} = identity_map{dim_ind} + d_source{dim_ind} .* options.integrator_b(stage_index);
            if ~isempty(options.integrator_b2)
                d_source2{dim_ind} = identity_map{dim_ind} + d_source2{dim_ind} .* options.integrator_b2(stage_index);
            end
            if (~options.single_sided)
                d_target{dim_ind} = identity_map{dim_ind} + d_target{dim_ind} .* options.integrator_b(stage_index);
                if ~isempty(options.integrator_b2)
                    d_target2{dim_ind} = identity_map{dim_ind} + d_target2{dim_ind} .* options.integrator_b2(stage_index);;
                end
            end
        end
        
        for dim_ind = 1:3
            b_source_deformation{dim_ind} = image_interp3(...
                b_source_deformation{dim_ind}, d_source ...
                , '*linear', 0, -1, window_mode);
            if (~options.single_sided)
                b_target_deformation{dim_ind} = image_interp3(...
                    b_target_deformation{dim_ind}, d_target ...
                    , '*linear', 0, -1, window_mode);
            end
            if ~isempty(options.integrator_b2)
                b2_source_deformation{dim_ind} = image_interp3(...
                    b2_source_deformation{dim_ind}, d_source2 ...
                    , '*linear', 0, -1, window_mode);
                if (~options.single_sided)
                    b2_target_deformation{dim_ind} = image_interp3(...
                        b2_target_deformation{dim_ind}, d_target2 ...
                        , '*linear', 0, -1, window_mode);
                end
            end
            % Prevent deformations from reaching outside the image:
            b_source_deformation{dim_ind} = max(...
                b_source_deformation{dim_ind}, 1);
            b_source_deformation{dim_ind} = min(...
                b_source_deformation{dim_ind}, source_size_xyz(dim_ind));
            if (~options.single_sided)
                b_target_deformation{dim_ind} = max(...
                    b_target_deformation{dim_ind}, 1);
                b_target_deformation{dim_ind} = min(...
                    b_target_deformation{dim_ind}, source_size_xyz(dim_ind));
            end
            if ~isempty(options.integrator_b2)
                b2_source_deformation{dim_ind} = max(...
                    b2_source_deformation{dim_ind}, 1);
                b2_source_deformation{dim_ind} = min(...
                    b2_source_deformation{dim_ind}, source_size_xyz(dim_ind));
                if (~options.single_sided)
                    b2_target_deformation{dim_ind} = max(...
                        b2_target_deformation{dim_ind}, 1);
                    b2_target_deformation{dim_ind} = min(...
                        b2_target_deformation{dim_ind}, source_size_xyz(dim_ind));
                end
            end
        end
        
        b_source_distance = b_source_distance + ...
            d_source_distance * options.integrator_b(stage_index);
        b_target_distance = b_target_distance + ...
            d_target_distance * options.integrator_b(stage_index);
        if ~isempty(options.integrator_b2)
            b2_source_distance = b2_source_distance + ...
                d_source_distance * options.integrator_b2(stage_index);
            b2_target_distance = b2_target_distance + ...
                d_target_distance * options.integrator_b2(stage_index);
        end
        
        clear d_source d_target d_source_distance d_target_distance
    end
    % Check if the step should be accepted and update the step size:
    maximum_step_size = max(options.time_span(2) - current_time, 0);
    step_size_controlled = false;
    any_too_large = 0;
    step_size_safety_factor = .8;
    % Two absolute error tolerances are implemented: distance and
    % deformation. Relative tolerance doesn't really make sense in
    % its usual form, especially during the first iteration where
    % only zero error is allowed.
    if (options.absolute_distance_tolerance > 0 && ~isempty(options.integrator_b2))
        step_size_controlled = true;
        b_source_distance_difference = ...
            abs(b2_source_distance - b_source_distance);
        b_target_distance_difference = ...
            abs(b2_target_distance - b_target_distance);
        
        any_too_large = any_too_large + ...
            b_source_distance_difference > ...
            options.absolute_distance_tolerance;
        any_too_large = any_too_large + ...
            b_target_distance_difference > ...
            options.absolute_distance_tolerance;
        
        distance_errors = [...
            b_source_distance_difference, ...
            b_target_distance_difference...
            ]
        % 0.8 is what some of Matlab's integration suite uses:
        maximum_step_size = min(...
            min(...
            current_step_size * ...
            (options.absolute_distance_tolerance ./ ...
            distance_errors).^(1./options.integrator_error_order) ...
            ) * step_size_safety_factor, maximum_step_size)
    end
    
    b_source_deformation_difference = [];
    b_target_deformation_difference = [];
    if ((options.absolute_deformation_tolerance > 0 && ~isempty(options.integrator_b2)))
        step_size_controlled = true;
        b_source_deformation_difference = b2_source_deformation;
        b_target_deformation_difference = b2_target_deformation;
        % This might be better done on velocities, or might it?
        for dim_ind = 1:3
            b_source_deformation_difference{dim_ind} = abs( ...
                b_source_deformation_difference{dim_ind} - ...
                b_source_deformation{dim_ind});
            if (~options.single_sided)
                b_target_deformation_difference{dim_ind} = abs( ...
                    b_target_deformation_difference{dim_ind} - ...
                    b_target_deformation{dim_ind});
            end
        end
    end
    
    
    if (options.absolute_deformation_tolerance > 0 && ~isempty(options.integrator_b2))
        step_size_controlled = true;
        % This might be better done on velocities, or might it?
        for dim_ind = 1:3
            any_too_large = any_too_large + sum(...
                b_source_deformation_difference{dim_ind} > ...
                options.absolute_deformation_tolerance);
            dd = options.absolute_deformation_tolerance ./ b_source_deformation_difference{dim_ind};
            maximum_step_size = min(...
                current_step_size * min(...
                dd .^ (1./options.integrator_error_order)) ...
                * step_size_safety_factor, maximum_step_size)
            %
            if (~options.single_sided)
                any_too_large = any_too_large + sum(...
                    b_target_deformation_difference{dim_ind} > ...
                    options.absolute_deformation_tolerance);
                dd = options.absolute_deformation_tolerance ./ ...
                    b_target_deformation_difference{dim_ind};
                maximum_step_size = min(...
                    current_step_size * min(...
                    dd .^ (1./options.integrator_error_order)) ...
                    * step_size_safety_factor, maximum_step_size)
            end
        end
    end
    
    % Control the deformation field's change using step size:
    if (any(options.maximum_deformation_per_step > 0) && ~isempty(options.integrator_b2))
        step_size_controlled = true;
        for dim_ind = 1:3
            if options.maximum_deformation_per_step(dim_ind) <= 0
                continue
            end
            source_deformation_ratio = ...
                options.maximum_deformation_per_step(dim_ind) / max(abs(b_source_deformation{dim_ind}(:) - current_source_deformation{dim_ind}(:)));
            any_too_large = any_too_large + ...
                (source_deformation_ratio < 1);
            maximum_step_size = min(...
                current_step_size * source_deformation_ratio ...
                * step_size_safety_factor, maximum_step_size);
            if (~options.single_sided)
                target_deformation_ratio = ...
                    options.maximum_deformation_per_step(dim_ind) / max(abs(...
                    b_target_deformation{dim_ind}(:) - current_target_deformation{dim_ind}(:)));
                any_too_large = any_too_large + ...
                    (target_deformation_ratio < 1);
                maximum_step_size = min(...
                    current_step_size * target_deformation_ratio ...
                    * step_size_safety_factor, maximum_step_size);
            end
        end
    end
    
    
    % Reduce step size if registration error goes up:
    if (options.registration_error_failure_step_scale > 0 && ...
            length(registration_error_array) > 1)
        step_size_controlled = true;
        registration_error_failures = ...
            registration_error_array(2:end) - ...
            registration_error_array(1:end - 1)...
            >= 0;
        if registration_error_failures(end)
            any_too_large = any_too_large + 1;
            maximum_step_size = min(...
                current_step_size * options.registration_error_failure_step_scale ...
                * step_size_safety_factor, maximum_step_size)
        end
    end
    
    % Step size remains constant if there are no error tolerances.
    next_step_size = current_step_size;
    if step_size_controlled
        next_step_size = maximum_step_size;
    end
    
    
    if (any_too_large > 0)
        % Step size is too large, so set it lower and recompute this
        % iteration:
        last_rejected_step_size = current_step_size;
        last_reject_iteration = iteration_index;
        number_rejected_steps = number_rejected_steps + 1;
        current_step_size = next_step_size;
        continue
    end
    
    % Update the deformations if the step is accepted:
    % Deformation deformation:
    
    current_source_deformation = b_source_deformation;
    current_target_deformation = b_target_deformation;
    current_source_distance = b_source_distance;
    current_target_distance = b_target_distance;
    
    current_source = image_interp3(...
        source, current_source_deformation ...
        , '*linear', 0, -1, window_mode);
    if (~options.single_sided)
        current_target = image_interp3(...
            target, current_target_deformation ...
            , '*linear', 0, -1, window_mode);
    end
    
    % Write images of the deformed source and target for quick
    % visual evaluation during these long computations:
    if options.save_intermediate_images
        deformed_source_image = current_source;
        deformed_source_image(isnan(deformed_source_image))=0;
        deformed_source_image = reshape(...
            contrast_stretch(deformed_source_image), ...
            size(deformed_source_image, 1), []);
        if (~options.single_sided)
            deformed_target_image = current_target;
            deformed_target_image(isnan(deformed_target_image))=0;
            deformed_target_image = reshape(...
                contrast_stretch(deformed_target_image), ...
                size(deformed_target_image, 1), []);
        end
        
        img = zeros(0, size(deformed_source_image, 2));
        if options.always_save_full_intermediates
            img = [img; reshape_contrast(source)];
            img = [img; reshape_contrast(target)];
        end
        img = [img; deformed_source_image];
        if (~options.single_sided)
            img = [img; deformed_target_image];
        end
        
        imwrite(img, ...
            [base_filename num2str(iteration_index, '_iter%05d') ...
            '_dfrms,dfrmt.png'] ...
            )
    end
    clear deformed_source_image deformed_target_image
    
    % Store iteration's state:
    if (options.save_intermediates)
        current_iteration_filename = [...
            base_filename ...
            num2str(iteration_index, '_iter%05d')...
            ];
        
        total_wall_time = toc(start_wall_time);
        total_cpu_time = cputime - start_cpu_time;
        if options.just_compute_distance
        else
            save([current_iteration_filename '.mat'], ...
                'current_source', 'current_target', ...
                'current_source_deformation', 'current_target_deformation', ...
                'current_source_distance', 'current_target_distance', 'number_rejected_steps', ...
                'total_cpu_time', 'total_wall_time' );
        end
        
        iteration_filenames = [iteration_filenames, ...
            {current_iteration_filename}];
    elseif ~options.just_compute_distance
        temp = struct();
        temp.current_source = current_source;
        temp.current_target = current_target;
        temp.current_source_deformation = current_source_deformation;
        temp.current_target_deformation = current_target_deformation;
        temp.current_source_distance = current_source_distance;
        temp.current_target_distance = current_target_distance;
        
        iteration_steps = [iteration_steps, temp];
    end
    
    
    % Update stuff:
    converged = converged || ...
        iteration_index >= options.maximum_iterations;
    iteration_index = iteration_index + 1;
    current_time = current_time + current_step_size;
    
    time_points = [time_points, current_time];
    source_distance_array = [source_distance_array, ...
        current_source_distance];
    target_distance_array = [target_distance_array, ...
        current_target_distance];
    
    
    % use known_distance and stop_source_distance to determine if no more interpolation is necessary:
    if options.known_distance >= 0 && options.stop_source_distance >= 0
        if options.verbose > 0
            fprintf('sd %f >= %f?\n', source_distance_array(end), options.stop_source_distance)
            fprintf('td %f >= %f?\n', target_distance_array(end), options.known_distance - options.stop_source_distance)
        end
        if source_distance_array(end) >= options.stop_source_distance || target_distance_array(end) >= options.known_distance - options.stop_source_distance
            break
        end
    end
    
    
    if options.convergence_tolerance > 0
        % Change in total distance added per unit time being less than
        % convergence_tolerance is considered convergence:
        source_distance_difference = source_distance_array(end) - ...
            source_distance_array(end - 1);
        target_distance_difference = 0;
        if (~options.single_sided)
            target_distance_difference = target_distance_array(end) - ...
                target_distance_array(end - 1);
        end
        convergence_ratio = (...
            (source_distance_difference + target_distance_difference) ./ ...
            (time_points(end) - time_points(end-1))) ./ (...
            source_distance_array(end) + target_distance_array(end));
        converged = converged || ...
            convergence_ratio <= options.convergence_tolerance;
    end
    
    registration_error = compute_registration_error(current_source, current_target);
    registration_error_array = [registration_error_array, ...
        registration_error];
    if (options.convergence_registration_error > 0)
        converged = converged || ...
            registration_error <= options.convergence_registration_error;
    end
    absolute_error = sum(abs(current_source(:)-current_target(:)));
    absolute_error_array = [absolute_error_array, ...
        absolute_error];
    if (options.convergence_absolute_error > 0)
        converged = converged || ...
            absolute_error <= options.convergence_absolute_error;
    end
    if (options.convergence_absolute_error_difference > 0 && length(absolute_error_array) > options.convergence_absolute_error_difference_window_size)
        converged = converged || ...
            -mean(diff(absolute_error_array(end - options.convergence_absolute_error_difference_window_size:end))) <= options.convergence_absolute_error_difference;
    end
    
    if (options.maximum_registration_error_failures > 0)
        number_registration_error_failures = sum(...
            registration_error_array(2:end) - ...
            registration_error_array(1:end - 1)...
            >= 0)
        excessive_registration_error_failures = ...
            number_registration_error_failures > ...
            options.maximum_registration_error_failures;
        converged = converged || excessive_registration_error_failures;
    end
    
    previous_step_size = current_step_size;
    current_step_size = next_step_size;
    if current_step_size <= 0
        break
    end
end


% Delete stages:
if ~options.just_compute_distance && options.save_stages
    for a = max(iteration_index - 2, 1):max(iteration_index - 1, 1)
        for prior_stage_index = 1:number_stages
            prior_stage_filename = [...
                base_filename ...
                num2str(a, '_iter%05d')...
                num2str(prior_stage_index, '_stage%05d')...
                ];
            delete([prior_stage_filename '.mat'])
        end
    end
end

r_structure.end_time = current_time;
r_structure.total_iterations = iteration_index - 1;
r_structure.source_deformation = current_source_deformation;
r_structure.source_distance = current_source_distance;
r_structure.source = current_source;
r_structure.total_wall_time = toc(start_wall_time);
r_structure.total_cpu_time = cputime - start_cpu_time;
r_structure.time_points = time_points;
r_structure.source_distance_array = source_distance_array;
r_structure.number_rejected_steps = number_rejected_steps;
r_structure.options = options;
r_structure.base_filename = base_filename;
r_structure.iteration_filenames = iteration_filenames;
r_structure.registration_error_array = registration_error_array;

if (~options.single_sided)
    r_structure.target_deformation = current_target_deformation;
    r_structure.target_distance = current_target_distance;
    r_structure.target = current_target;
    r_structure.target_distance_array = target_distance_array;
end

if ~options.just_compute_distance && ~options.save_intermediates
    r_structure.iteration_steps = iteration_steps;
end

warning(prior_warning_state);
end