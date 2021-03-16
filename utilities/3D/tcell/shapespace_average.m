function [I_out, r_out] = shapespace_average(input_images, input_image_names, t_cell_info, dist_param)
%
% the method is based on the computing of diffeomorphic distances and do
% embedding to get the coordinates of the cells in the shapespace and then
% compute the average position. At last, resample a shape based on the
% coordinate. 
%
%
% Author: Xiongtao Ruan 
% Date: Oct. 30, 2015


% Set random seed:
    

s = RandStream('mt19937ar', 'Seed', 892734);
RandStream.setGlobalStream(s);

registration_options = dist_param.registration_options;
working_dir = dist_param.working_dir;
compute_columns = dist_param.compute_columns;

average_filename = [working_dir, 'average_shape_info.mat'];

if exist(average_filename, 'file')
    load (average_filename);
    return;
end

dist_dir = [working_dir, 'distances/'];
mkdir_recursive(dist_dir);

options = struct();
options.embedding = struct();
options.embedding = process_options_structure( struct( ...
                            'method', 'most_complete', ...
                            'force_positive_definiteness', true, ... %% this is only for the embed_partial_distance_matrix function
                            'explained_variance_threshold', 0.9, ...
                            'maximum_shape_space_dimension', 50 ...
                            ), options.embedding);
                        
embedding_methods = {'regress', 'regress landmark', 'mdscale'};
options.embedding.method = embedding_methods{1};


condition_retrieval_function = t_cell_info.synapse_info.condition_retrieval_function;
sensor_retrieval_function = t_cell_info.synapse_info.sensor_retrieval_function;
run_date_retrieval_function = t_cell_info.synapse_info.run_date_retrieval_function;
run_num_retrieval_function = @(given_strings)cellfun(@(x)x{1}{1},regexp(given_strings,'[^_/]+ - [^_/]+_[^_/]+_[^_/]*[0-9]{2} *[0-9]{2} *[0-9]{2} *[^_/]* *[Rr]un *([0-9]+) *\.csv','tokens'),'UniformOutput',false);
frame_num_retrieval_function = @(given_strings)cellfun(@(x)x{1}{1},regexp(given_strings,'[^_/]+ - [^_/]+_[^_/]+_[^_/]+_frame000([0-9]{2})_[^_/]+','tokens'),'UniformOutput',false);
synapse_retrieval_function = @(given_strings)cellfun(@(x)[x{1}{1}, x{1}{2}], regexp(given_strings, '[^_/]+ - [^_/]+_[^_/]+_[^_/]+_frame000[0-9]{2}_synapse*00([0-9]{3}), *00([0-9]{3})*','tokens'),'UniformOutput',false);                                                                                                 
img_num = numel(input_images);

conditions = condition_retrieval_function(input_image_names);
sensors = sensor_retrieval_function(input_image_names);
run_dates = run_date_retrieval_function(input_image_names);
run_nums =  run_num_retrieval_function(input_image_names);
frame_nums = frame_num_retrieval_function(input_image_names);
synapses = synapse_retrieval_function(input_image_names);

short_condition_set = {'F', 'B', 'R'};
condition_types = {'Full Stimulus', 'Blockade', 'active Rac Cofilin'};
short_conditions = cellfun(@(x) short_condition_set(~cellfun(@isempty, cellfun(@(y) regexp(x, y), condition_types, 'UniformOutput', false))), conditions);
short_names = cellfun(@(sc, ss, srn, srd, sfn, ssy) sprintf('%s_%s_R%s_%d%02d%02d_f%s_%s', sc, ss, srn, srd(1), srd(2), srd(3), sfn, ssy), short_conditions, sensors, run_nums, run_dates, frame_nums, synapses, 'UniformOutput', false);

[short_names, inds] = sort(short_names);
input_image_names = input_image_names(inds);
input_images = input_images(inds);

rand_order = randperm(img_num);
used_columns = rand_order(1 : compute_columns);


compute_ind_combine = zeros(img_num) - diag(ones(img_num, 1));
compute_ind_combine(:, used_columns) = compute_ind_combine(:, used_columns) + 1;
compute_ind_combine = compute_ind_combine > 0;
N_compute = sum(compute_ind_combine(:));
[y, x] = ind2sub([img_num, img_num], find(compute_ind_combine == 1));

dist_mat = -ones(img_num) + eye(img_num);
dist_mat_filename = [working_dir, 'dist_mat.mat'];
if exist(dist_mat_filename, 'file')
    load(dist_mat_filename);
end


computing_incomplete = true;
compute_average = true;

dir_info = dir_recursive(dist_dir);

if all(cellfun(@(x) sum(~cellfun(@isempty, regexp(dir_info, x))) == img_num - 1, short_names(used_columns)))
    computing_incomplete = false;
end


if compute_average 
    % compute distances
    while computing_incomplete
        r = struct();
        for i = 1 : N_compute
            x_i = x(i);
            y_i = y(i);
            if x_i < y_i
                s_i = x_i;
                t_i = y_i;
            else
                s_i = y_i;
                t_i = x_i;
            end
            source_short_name = short_names{s_i};
            target_short_name = short_names{t_i};

            dir_i = [dist_dir, source_short_name, '/'];
            if ~exist(dir_i, 'dir')
                mkdir(dir_i);
            end

            temp_filename = sprintf('%s%s-%s.tmp', dir_i, source_short_name, target_short_name);
            dist_filename = sprintf('%s%s-%s.mat', dir_i, source_short_name, target_short_name);

            if exist(dist_filename, 'file') || exist(temp_filename,'file')
                % disp('The computing is finished or being working somewhere else');
                continue;
            else    
                fclose(fopen(temp_filename, 'w'));
            end        
            disp(sprintf('compute distance between cell %d and cell %d', s_i, t_i));
            source = input_images{s_i};
            target = input_images{t_i};

            % r = Greedy3D_lambda_pre_compressed(source, target, 1, registration_options)
            r_raw = Greedy3D_lambda_pre(source, target, 1, registration_options)
            
            r.total_distance = r_raw.total_distance;
            save(dist_filename, 'registration_options', 'r');
            if exist(dist_filename, 'file') && exist(temp_filename, 'file')
                delete(temp_filename);
            end
        end

        % delete old temp files

        dir_info = dir_recursive(dist_dir);

        is_temp_file = ~cellfun(@isempty, regexp(dir_info, '.tmp')); 
        if sum(is_temp_file) > 0 
            temp_filenames = dir_info(is_temp_file);
            for i = 1 : numel(temp_filenames)
                current_temp_filename = temp_filenames{i};
                temp_file_info = dir(current_temp_filename);
                if (datenum(datetime('now')) - [temp_file_info.datenum]) * 24 * 60 > 10
                    delete(current_temp_filename);
                end
            end
        end

        if all(cellfun(@(x) sum(~cellfun(@isempty, regexp(dir_info, x))) == img_num - 1, short_names(used_columns)))
            computing_incomplete = false;
            continue;
        end    
        pause(30);
    end
    
    post_dist_process_flag = [working_dir, 'post_dist_processing.tmp'];
    
    if ~exist(post_dist_process_flag, 'file') && ~exist(average_filename, 'file')
        fclose(fopen(post_dist_process_flag, 'w'));
        
        dist_mat_filename = sprintf('%sdistance_matrix.mat', working_dir);
        
        if exist(dist_mat_filename, 'file')
            load(dist_mat_filename);
        end       

        %collecting all distances. 
        collecting_distances = true;
        if collecting_distances 
            N_save = 0;
            for i = 1 : N_compute
                x_i = x(i);
                y_i = y(i);
                if dist_mat(x_i, y_i) > 0 
                    continue;
                end
                if x_i < y_i
                    s_i = x_i;
                    t_i = y_i;
                else
                    s_i = y_i;
                    t_i = x_i;
                end

                source_short_name = short_names{s_i};
                target_short_name = short_names{t_i};

                dist_filename = sprintf('%s%s/%s-%s.mat', dist_dir, source_short_name, source_short_name, target_short_name)

                try
                    load(dist_filename);
                catch ME
                    if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
                        disp('cannot read the file, and delete the file');
                        delete(dist_filename);
                        [r] = compute_single_missing_distance(short_names, input_images, dist_filename, registration_options);
                    end
                end

                dist_mat(x_i, y_i) = r.total_distance;
                dist_mat(y_i, x_i) = dist_mat(x_i, y_i);
                N_save = N_save + 1;

                if rem(N_save, 10000) == 0 || i == N_compute
                    save(dist_mat_filename, 'dist_mat');
                end
            end
        end
        embedding_model_filename = sprintf('%sembedding_model.mat', working_dir);
        
        if exist(embedding_model_filename, 'file')
            load(embedding_model_filename);
        else
            dist_mat(dist_mat == -1) = NaN;
        % get the embedded positions
        
            [ embedding_model ] = embed_distance_matrix( dist_mat, options );
            save(embedding_model_filename, 'embedding_model');
            clear dist_mat;
        end
        positions = embedding_model.positions;
        
        mean_position = mean(positions);
        sample_average_shape = true;
        
        % sample the average image
        if sample_average_shape
            image_function = @(x) input_images{x};
            shape_space2 = {embedding_model.positions, embedding_model.convex_hull, embedding_model.tessellation};
            r_out = render_point_3d_windowed(mean_position, shape_space2, image_function, 0, registration_options); 
            I_out = r_out.interpolated_image.image{1};
            I_out = contrast_stretch(round(I_out));
            save(average_filename, 'r_out', 'positions', 'I_out');
            delete(post_dist_process_flag);
        end
        
    end
                  
    while exist(post_dist_process_flag, 'file')
        pause(30);
    end
    load(average_filename);
end
     
 
        
    

end


function [r] = compute_single_missing_distance(short_names, input_images, missing_filename, registration_options)
    %this function is used for a missing distance where the distance file
    %is corrupted. 
    img_inds = find(~cellfun(@isempty, cellfun(@(x) regexp(missing_filename, x), short_names, 'UniformOutput', false)) == 1);
    x_i = img_inds(1);
    y_i = img_inds(2);
    if x_i < y_i
        s_i = x_i;
        t_i = y_i;
    else
        s_i = y_i;
        t_i = x_i;
    end
    source_short_name = short_names{s_i};
    target_short_name = short_names{t_i};

    source = input_images{s_i};
    target = input_images{t_i};

    % r = Greedy3D_lambda_pre_compressed(source, target, 1, registration_options)
    r = Greedy3D_lambda_pre(source, target, 1, registration_options)

    save(missing_filename, 'registration_options', 'r');
            

end







