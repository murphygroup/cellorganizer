function frame = makeFrame(model,options,curr_point,framefolder,currframe)
%Moved from model2diffeomorphicInstance.m
%This code generates a single frame from a diffeomorphic shape space model.
%Inputs: 
% model - a cellorganizer diffeomorphic model
% options - the set of options for synthesizing from a trained shape space
% model
% curr_point - the location in the shape space to synthesize
% framefolder - the location to save the resulting image
% currframe - the index of the current frame for book keeping.
%
%Outputs:
%Saves .tif files for cell and nuclear shape in the framefolder
%Save a .mat file containing the cells willmore energy
%

%Created by:
%Taraz Buck
%Moved to its own function by:
%Devin Sullivan 7/22/13
%
%Edited by: 
%D. Sullivan 9/6/13 - supressed warnings for willmore energy calculation
%D. Sullivan 9/12/13 - added check if result files were present to skip
%interpolation. 

%icaoberg 11/1/2012
try
   default_options.method = param.synthesis.diffeomorphic.method;
catch
   %icaoberg
   %options.method = 'convnfft';
   default_options.method = 'convnfft_fast';
end

%yy 11/28/2012
try
   default_options.convergence_absolute_error = param.synthesis.diffeomorphic.convergence_absolute_error;
catch

end

%icaoberg 11/9/2012
%maximum number of iterations
try
   default_options.maximum_iterations = param.synthesis.diffeomorphic.maximum_iterations;
catch
   default_options.maximum_iterations = 100;
end

%changes per step
try
   default_options.convergence_tolerance = param.synthesis.diffeomorphic.convergence_tolerance;
catch
   default_options.convergence_tolerance  = 0;
   %options.convergence_tolerance  = 8.279212e-02;
end

%changes per step
try
   default_options.convergence_registration_error_scale = param.synthesis.diffeomorphic.convergence_registration_error_scale;
catch
   default_options.convergence_registration_error_scale  = 5e-3;
end

%tebuck 11/5/2012
%forward euler method as rk integrator
%options.integrator_c = 0;
%options.integrator_b = 1;
%options.integrator_a = 0;
%options.integrator_b2 = [];
%options.integrator_first_same_as_last = false;    

%tebuck 11/9/2012
default_options.alpha = 1;          % coefficient of Lapacian operator 
default_options.gamma = 0.04;       % Coefficient of the Identity 
default_options.keep_proportion = 0; 
default_options.single_sided = false; 
default_options.max_time = 50000; 
default_options.step_size = 0.01; 
default_options.use_fft = true; 
default_options.use_gaussian_kernel = false; 
default_options.absolute_distance_tolerance = 0; 
default_options.absolute_deformation_tolerance = 0; 
default_options.periodic_space = false; 
default_options.use_compression = false; 
default_options.drop_kernel = true; 
default_options.maximum_registration_error_failures = 0;

%tebuck 11/25/2012
%note that the following code is intended for use with values of 1 or 2
%downsampling scale
try
  default_options.downsampling_scale = param.synthesis.diffeomorphic.downsampling_scale;
catch
  default_options.downsampling_scale = 1;
end

default_options.maximum_deformation_per_step = [4/default_options.downsampling_scale, 4/default_options.downsampling_scale, 0.5];

%tebuck updated radius to match new model
default_options.window_radius = 192/default_options.downsampling_scale; 
default_options.filter_radius = 16/default_options.downsampling_scale;

options = process_options_structure(default_options, options);


celldir = [framefolder filesep 'Cellwalk'];
if ~isdir(celldir)
    mkdir(celldir);
end
nucdir = [framefolder filesep 'Nucwalk'];
if ~isdir(nucdir)
    mkdir(nucdir);
end
edir = [framefolder filesep 'energies'];
if ~isdir(edir)
    mkdir(edir);
end

cellframe = [celldir filesep 'frame' num2str(currframe) '.tif'];
nucframe = [nucdir filesep 'frame' num2str(currframe) '.tif'];
eframe = [edir filesep 'walkinfo' num2str(currframe) '.mat'];


if exist(cellframe,'file')&&exist(nucframe,'file')
    disp('Cell and Nuclear interpolations found. Skipping recalculation.')
    if exist(eframe,'file')
        return
    else
        energy = AverageEnergy(cellimg);
        save(eframe,'energy');
        return
    end
end
    
options.tempparent = framefolder;


framefile =  [framefolder filesep 'frame_info.mat'];

if ~exist(framefile, 'file')
    frame = generate_frame_from_trained_shape_space_model(model, curr_point, options );
    save(framefile, 'frame')
else
    disp(['Frame file ' framefile ' already exists. Loading.']);
    load(framefile)
end

%icaoberg 10/8/2012
%imwrite(reshape_contrast(round(frame.interpolated_image.get_data())), image_filename)
img = reshape_contrast(round(frame.interpolated_image.get_data()));
img = reshape( img, size(img,1), size(img,1), [] );

%tebuck 11/5/2012
%the variable image must contain three unique values: 0 for background, 1 for cell and
%2 for nucleus, i.e. img is an indexed image

values = unique( img );
nucimg = zeros( size(img) );
%D. Sullivan 12/15/14 - need to correct for when nuc is not the pre-assumed
%level or not present! 
if isfield(model.nuclearShapeModel,'level')
    %If it's 0 or nan, assume the nuclear model does not exist. 
    if model.nuclearShapeModel.level~=0 && ~isnan(model.nuclearShapeModel.level)
        %(add 1 for matlab, non 0 indexing)
        nucimg(find(img == values(model.nuclearShapeModel.level+1))) = 255;
    end
else
    warning('No level found for nuclear model, assuming nucleus is level 3')
    nucimg(find(img == values(end))) = 255;
end
    


cellimg = zeros( size(img) );
%icaoberg 11/5/2012
%D. Sullivan 12/15/14 - need to correct for when the cell is not the
%pre-assumed level or not present! 
if isfield(model.nuclearShapeModel,'level')
    %If it's 0 or nan, assume the nuclear model does not exist. 
    if model.cellShapeModel.level~=0 || ~isnan(model.cellShapeModel.level)
        %(add 1 for matlab, non 0 indexing)
        cellimg(find(img == values(model.cellShapeModel.level+1))) = 255;
    end
else
    warning('No level found for nuclear model, assuming cell is level 2')
    cellimg(find(img == values(2))) = 255;
end

%D. Sullivan 6/27/13 - added Willmore energy tracking so we can
%directly observe the cell energy for the current frame
%         energy(i) = AverageEnergy(frame{i}.interpolated_image.image{1});
%D. Sullivan 9/6/13 suppressed warnings for this step
warning('off','all')
energy = AverageEnergy(cellimg);
warning('on','all')

%D. Sullivan 2/25/13 need to remove the padding zeros from the
%framework model maintaining a 1 slice padding (this is for converting
%to meshing if needed later.
keepslices = find(sum(sum(cellimg,1),2)~=0);

cellimg_temp = cellimg; 
nucimg_temp = nucimg;
if min(keepslices) == 1
    cellimg_temp = padarray(cellimg_temp, [0,0,1], 'pre');
    nucimg_temp = padarray(nucimg_temp, [0,0,1], 'pre');
end

if max(keepslices) == size(cellimg,3)
    cellimg_temp = padarray(cellimg_temp, [0,0,1], 'post');
    nucimg_temp = padarray(nucimg_temp, [0,0,1], 'post');
end

cellimg = cellimg_temp;
nucimg = nucimg_temp;

keepslices = find(sum(sum(cellimg,1),2)~=0);

%D. Sullivan 6/25/13 creation of cell array to save each cell/nuc
%in a random walk.
%         cellimg = cellimg(:,:,keepslices(1)-1:keepslices(end)+1);
%         nucimg = nucimg(:,:,keepslices(1)-1:keepslices(end)+1);

% is this pushed
%D. Sullivan 6/25/13 temporarily create a cellimg and nucimg tiff
%in a temp result folder for the random walk stuff
%D. Sullivan 7/16/ making tempfolder a part of param struct and
%renaming framefolder
%         tempfolder = './temp';

output_folder = [framefolder filesep 'Cellwalk'];
if ~exist( output_folder )
    mkdir( output_folder );
end
%Add check to allow for synthesis of 2D images
if ndims(cellimg)==3
    cellimg = cellimg(:,:,keepslices(1)-1:keepslices(end)+1);
end
img2tif(cellimg,[framefolder filesep 'Cellwalk' filesep 'frame' num2str(currframe) '.tif'])

output_folder = [framefolder filesep 'Nucwalk'];
if ~exist( output_folder )
	mkdir( output_folder )
end
%Add check to allow for synthesis of 2D images
if ndims(nucimg)==3
    nucimg = nucimg(:,:,keepslices(1)-1:keepslices(end)+1);
end
img2tif(nucimg,[framefolder filesep 'Nucwalk' filesep 'frame' num2str(currframe) '.tif'])

output_folder = [framefolder filesep 'energies'];
if ~exist( output_folder )
	mkdir( output_folder )
end

save([framefolder filesep 'energies' filesep 'walkinfo'...
    num2str(currframe) '.mat'],'energy');


