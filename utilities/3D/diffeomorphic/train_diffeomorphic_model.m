function [ shape_space_model ] = train_diffeomorphic_model(options_in)
%Main wrapper for Taraz's LDDMM shapespace training function
%default parameters take from test_shape_space_construction_synthesis
%

% Aug 30, 2013 G. Johnson    Changed they way files are input into the
%                            diffeomorphic model function
% Nov 4, 2013 G. Johnson     Set aspect ratio to [1,1,1] by default

    
    options_in = ml_initparam(options_in, struct( ...
                                            'imfunc', @(x) [options_in.preprocessingFolder filesep 'cell' num2str(x) '.mat'], ...
                                            'model', struct()));
                                        
    options_in.model = ml_initparam(options_in.model, struct('diffeomorphic', struct()));
    
    options_in.model.diffeomorphic = ml_initparam(options_in.model.diffeomorphic, ...
        struct( 'use_distance_matrix_completion', true, ...
                'minimum_relative_semidiameter', 1 / 4, ...
                'maximum_relative_semidiameter', 2 / 3, ...
                'generate_cycle', true, ...
                'useCurrentResults', false, ...
                'tempdir', [options_in.tempparent filesep 'diffeomorphic'], ...
                'downsample', max(options_in.model.resolution)./options_in.model.resolution, ...
                'com_align', true, ...
                'z_align', 'com', ... %can be 'bottom', 'com' or 'largest_slice'
                'number_windows', 1, ...
                'regfunc', @(x,y) cellfile2registeredimg(options_in.imfunc(x), y), ...
                'relative_filter_radius', 0.10, ...
                'scale_z_kernel', false ... %flag to scale the z kernel to be imsize(1)/imsize(3);
                ));
    
    options_in.model.diffeomorphic.resolution = options_in.model.diffeomorphic.downsample .* options_in.model.resolution;
        


            
    if ~exist(options_in.model.diffeomorphic.tempdir, 'dir')
        mkdir(options_in.model.diffeomorphic.tempdir)
    end
            
    options = struct();
    
    options.tempparent = options_in.tempparent;
    options.save_location = [options_in.model.diffeomorphic.tempdir];
    
    [imfunc_return_aligned_images, min_im_size, maxdims, image_output_size] = get_diffeo_image_function(options_in);
        
    image_width = min(min_im_size(1:2));
    
    options.number_images = options_in.documentation.numimgs;
     
    options.voxel_size = options_in.model.resolution;
    options.shape_aspect_ratio = [1, 1, 1];   
    
    %use filter that is 1/10 the size of the image
    options.filter_radius = ceil(image_output_size(1)*options_in.model.diffeomorphic.relative_filter_radius);
    
    %if the window radius is not set, we just use one window
    options.window_radius = ceil(image_output_size(1)/options_in.model.diffeomorphic.number_windows);
%   options.window_radius = max(min_im_size(1:2)); %old code

    
    %use z_radius kernel that is 1/5 the size of the z dimension
    if(maxdims(3) > 1)
        
        if options_in.model.diffeomorphic.scale_z_kernel
            options.kernel_z_radius = options.filter_radius;
        else
            options.kernel_z_radius = ceil(min_im_size(3)/5);
        end
    else
        options.kernel_z_radius = 0;
    end
 
    options.convergence_absolute_error = (mean([options_in.model.diffeomorphic.minimum_relative_semidiameter, options_in.model.diffeomorphic.maximum_relative_semidiameter]) * image_width * 2 * pi) / 10;
    options.maximum_deformation_per_step = ones(1, 3) .* [1,1,0.5]; 
%     options.maximum_deformation_per_step = ones(1, 3) * .25; 
    options.convergence_absolute_error_difference = 25;
    options.useCurrentResults = options_in.model.diffeomorphic.useCurrentResults;
    
    if options_in.model.diffeomorphic.use_distance_matrix_completion
        options.desired_shape_space_dimensionality = 2;
        % shape_space_options.desired_shape_space_dimensionality = 5;
    end
    
    %Resize the images to be divisible by the window_raidus
    options.image_function = @(x) pad_imfunc_to_window_size(x, imfunc_return_aligned_images, options.window_radius);
    options.imsize = min_im_size;
    
    
    options = ml_initparam(options_in.model.diffeomorphic, options);
    
    shape_space_model = train_shape_space_model(options);

    
%     shape_space_model.imfunc = imfunc_return_aligned_images;
%     shape_space_model.numimgs = param.documentation.numimgs;
%     shape_space_model.imsize = min_im_size;
%     shape_space_model.name = 'diffeomorphic model of the cell';
%     shape_space_model.type = 'diffeomorphic';
% %     shape_space_model.matCompletionFunctionString = options.matCompletionFunctionString;
%     shape_space_model.version = 1.0;

end


