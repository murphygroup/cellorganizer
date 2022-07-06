function answer = img2slml( varargin )
% IMG2SLML Trains a generative model of subcellular location from a
% collection of images and saves the model to disk.
%
% A CellOrganizer model consists of four components,
%
% 1) a (optional) documentation component
% 2) a nuclear membrane model,
% 3) a cell membrane model and,
% 4) a protein pattern model.
%
% ┌────────────────────────┐
% │List Of Input Parameters│
% └────────────────────────┘
% Inputs                      Descriptions
% ------                      ------------
% dimensionality              2D/3D
% dnaImagesDirectoryPath      DNA images collection directory, list of files or pattern
% cellImagesDirectoryPath     Cell images collection directory, list of files or pattern
% proteinImagesDirectoryPath  Protein images collection directory, list of files or pattern
% options                     List of options
%
% The input argument options holds the valid parameters for all of these components.
%
% ┌───────────────┐
% │List Of Options│
% └───────────────┘
% Mandatory options         Descriptions
% -----------------         ------------
% model.resolutions         Any double 1x2/1x3 double vector.
%                           (microns/voxel).
%
% Generic model options     Descriptions
% ---------------------     ------------
% masks                     (optional) Masks collection directory.
%
% train.flag                (optional) Selects what model is going to be trained ('nuclear',
%                           'framework', or 'all'). Default is 'all'.
%
% model.name                (optional) Holds the name of the model. Default is empty.
% model.id                  (optional) Holds the id of the model. Default is a randomly generated string.
% model.filename            Holds the output filename.
% downsampling              Downsampling vector to be used during preprocessing.
%
% Debugging options
% -----------------
% debug                     If set to true, then the function will (1) keep temporary results folder, (2) will
%                           print information useful for debugging. Default is false.
% display                   If set to true, then plots useful for debugging with be open. This functionality is
%                           meant for debugging only, setting this to true will considerably slow down
%                           computation. Default is false;
% save_segmentations        Will save the segmentations to the model file. Setting this option to true will create
%                           a considerably large file.
%
% Nuclear shape model options  Descriptions
% ---------------------------  ------------
% nucleus.class                (mandatory) Holds the nuclear membrane model class.
% nucleus.type                 (mandatory) Holds the nuclear membrane model type.
% nucleus.name                 (optional) Holds the name of the nuclear model. Default is empty.
% nucleus.id                   (optional) Holds the id of the nuclear model. Default is a randomly generated string.
%
% Cell shape model options  Descriptions
% ------------------------  ------------
% cell.class                (mandatory) Holds the cell membrane model class.
% cell.type                 (mandatory) Holds the cell membrane model type.
% cell.name                 (optional) Holds the name of the cell model. Default is empty.
% cell.id                   (optional) Holds the id the cell model. Default is empty.
%
% Protein shape model options  Descriptions
% ---------------------------  ------------
% protein.class                (mandatory) Holds the protein model class.
% protein.type                 (mandatory) Holds the protein model type.
% protein.name                 (optional) Holds the name of the protein model. The default is empty.
% protein.id                   (optional) Holds the id of the protein model. Default is a randomly generated string.
% protein.cytonuclearflag      (optional) Determines whether the protein pattern will be generated in
%                              the cytosolic space ('cyto'), nuclear space ('nuc') or everywhere ('all').
%                              Default is cyto.
%
% ┌────────────────────────────────────┐
% │List Of Options per Model class/type│
% └────────────────────────────────────┘
% 2D PCA model options
% --------------------
% model.pca.latent_dim      (optional) This specifies how many latent dimensions should be used for modeling
%                           the shape space. Valid values arepositive integers. The default is 15.
%
% 2D diffeomorphic model options
% ------------------------------
% model.diffeomorphic.distance_computing_method     (optional) ‘faster'
% model.diffeomorphic.com_align                     (optional) 'nuc'
%
% T cell distribution model options
% ---------------------------------
% model.tcell.synapse_location            (mandatory) File path to annotation of the synapse positions of the T cells as input.
% model.tcell.results_location            (mandatory) File path for where the results should be saved.
% model.tcell.named_option_set            (mandatory) The running choice for CellOrganizer and one sensor of two-point annotation.
% model.tcell.use_two_point_synapses      (optional) Set up the mode of synapse to use, as a default, we use one-point,
%                                         if needed you can use two-point by set up the option as true.
% model.tcell.sensor                      Set up protein name.
% model.tcell.timepoints_to_include       (optional) If creation of models for only a subset of the time points is desired,
%                                         edit to specify which time points to include.
% model.tcell.model_type_to_include       (mandatory) Set up for model to include.
% model.tcell.infer_synapses              (mandatory) true or false.
% model.tcell.adjust_one_point_alignment  (optional) Set up alignment adjustment true or false.
% model.tcell.ometiff                     (optional) If true, then it assumes images are OME.TIFFs with annotations. Default is false.
%
% 3D SPHARM-RPDM model options
% ----------------------------
% model.spharm_rpdm.alignment_method 	(optional) method by which cells willbe aligned when producing shape descriptors
%                                       The possible values are 'major_axis' (defaut) or 'foe'.
% model.spharm_rpdm.rot ation_plane 	(optional)  Dimensions of image that will used for alignment.
%                                       The possible values are 'xy' (defaut), 'xz', 'yz' or ‘xyz'. For example,
%                                       if ‘xy‘ is specified, each cellwill be rotated in the 	xy plane (around the z axis).
% model.spharm_rpdm.postprocess         (optional) This specifies whether alignment and size normalization
%                                       should be done after parameterization.  The values are ‘true’ (default) and ‘false’.
% model.spharm_rpdm.maxDeg              (optional) This specifies the degree up to which spherical harmonics
%                                       should be calculated.  Valid values are positive integers.  The default is 31.
% model.spharm_rpdm.components          (mandatory) This specifies which components should be included in the
%                                       shape model.  The valid values are {'cell'}, {'nuc'}, or {'cell', 'nuc'}.
% model.spharm_rpdm.latent_dim          (optional) This specifies how many latent dimensions should be used for
%                                       modeling the shape space.  Valid values are positive integers.  The default is 15.
%
% ┌─────────────┐
% │Documentation│
% └─────────────┘
% This is an optional structure with multiple elements that holds documentation about this model.
%
% documentation.<name>      Holds the value of variable <name>. This is meant to be meta information. Default is empty.
%
% Helper Options
% -------------
% verbose                   (optional) Displays messages to screen. Default is true.
% debug                     (optional) Reports errors and warnings. Default is false.

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2007-2020 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://www.cellorganizer.org or
% send email to cellorganizer@compbio.cmu.edu

answer = false;

if ~exist( [ pwd filesep 'log'], 'dir' )
    mkdir( 'log' );
end

logfile = datestr(datetime('now', 'TimeZone', 'local'), 'yyyymmddHHMMSS.FFF'); %#ok<AGROW>
logfile = [ pwd filesep 'log' filesep logfile, '.log' ];
diary( logfile )

if isdeployed
    disp('Running deployed version of img2slml');
    
    disp('Checking number of input arguments')
    if length(varargin) == 1
        text_file = varargin{1};
    else
        error('Deployed function takes only 1 argument. Exiting method.');
        return
    end
    
    disp('Checking existence of input file')
    [filepath, name, ext] = fileparts(text_file);
    
    if ~exist(text_file, 'file')
        warning('Input file does not exist. Exiting method.');
        return
    end
    
    disp(['Attempting to read input file ' text_file]);
    fid = fopen(text_file, 'r' );
    
    disp('Evaluating lines from input file');
    while ~feof(fid)
        line = fgets(fid);
        disp(line);
        try
            eval(line);
        catch err
            disp('Unable to parse line');
            getReport(err)
            return
        end
    end
    
    disp('Closing input file')
    fclose(fid);
    
    disp( 'Checking nuclear membrane images list' )
    if ~isempty( dnaImagesDirectoryPath )
        dnaImagesDirectoryPath = parse_ometiff_deployed(dnaImagesDirectoryPath);
    else
        disp( 'List is empty' )
    end
    
    disp( 'Checking cell membrane images list' )
    if ~isempty( cellImagesDirectoryPath )
        cellImagesDirectoryPath = parse_ometiff_deployed(cellImagesDirectoryPath);
    else
        disp( 'List is empty' )
    end
    
    disp( 'Checking protein pattern images list' )
    if ~isempty( proteinImagesDirectoryPath )
        proteinImagesDirectoryPath = parse_ometiff_deployed(proteinImagesDirectoryPath);
    else
        disp( 'List is empty' )
    end
    
    disp( 'Checking mask pattern images list' )
    if ~isempty( options.masks )
        options.masks = parse_ometiff_deployed(options.masks);
    else
        disp( 'List is empty' )
    end
    
    disp('Checking and getting default parameters')
    param = get_cellorganizer_default_parameters( 'training', options );
else
    disp('Running img2slml')
    disp('Check number of input arguments')
    try
        param = varargin{5};
    catch
        param = [];
    end
    
    dimensionality = varargin{1};
    dnaImagesDirectoryPath = varargin{2};
    cellImagesDirectoryPath = varargin{3};
    proteinImagesDirectoryPath = varargin{4};
    param = varargin{5};
    
    disp('Checking and getting default parameters')
    param = get_cellorganizer_default_parameters( 'training', param );
end

%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK INPUT ARGUMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' '); print_large_title('Validating input arguments' );
if isdeployed || length(varargin) == 5
    [dnaImagesDirectoryPath, cellImagesDirectoryPath, ...
        proteinImagesDirectoryPath,labels] = check_images_directory_paths( ...
        dnaImagesDirectoryPath, cellImagesDirectoryPath, ...
        proteinImagesDirectoryPath, param );
    
    if ~ismember( param.train.flag, {'protein'} ) && ...
            isempty( dnaImagesDirectoryPath ) && ...
            isempty( cellImagesDirectoryPath ) && ...
            isempty( proteinImagesDirectoryPath )
        answer = false;
        return
    end
    
    if (isfield( param, 'protein' ) && ...
            isfield( param.protein, 'class' ) && ...
            strcmpi( param.protein.class, 'standardized_voxels' )) && ...
            (isfield( param, 'protein' ) && ...
            isfield( param.protein, 'type' ) && ...
            strcmpi( param.protein.type, 'standardized_map_half-ellipsoid' ))
        disp(' '); print_large_title('Cleaning up data' );
        if isfield(param.model, 'tcell') && isfield(param.model.tcell, 'ometiff') && ...
                (param.model.tcell.ometiff == true)
            get_tif_annotation_from_ometiff( proteinImagesDirectoryPath );
        end
    end
    param = clean_up_training_input_arguments( dimensionality, param );
    if isempty( param )
        warning('Unable to extract one or several parameters. Exiting method');
        return
    end
else
    error('Wrong number of input arguments. Exiting method.');
    return
end

if param.verbose
    fprintf( 1, '%s\n', 'Checking the existence of temporary folder' );
end

% %TRAIN THE GENERATIVE MODEL
disp(' '); print_large_title('Training generative model' );
param.dimensionality = dimensionality;
model = img2model( dnaImagesDirectoryPath, ...
    cellImagesDirectoryPath, ...
    proteinImagesDirectoryPath, ...
    param );

if isempty( model ) || numel(fields(model)) <=2
    warning( ['Method img2model returned an empty model. Exiting method.']);
    answer = false;
    return
end

%PARSE GENERATIVE MODEL INTO SLML INSTANCE
disp(' '); print_large_title('Parse and clean generative model' );
if strcmpi( dimensionality, '2D' )
    disp(upper('Adding parameters to model structure'));
    model = parse_and_clean_generative_model( model, param );
    
    disp(upper('Adding documentation to model structure'));
    model = add_documentation_to_model( model, param );
    
    disp(upper('Clean up and wrap up model structure'));
    model = clean_up_and_wrap_up_model( model, param );
    
    disp(' '); disp(upper('Clean up workspace and environment'));
    answer = clean_up( param );
else
    disp(upper('Adding parameters to model structure'));
    model = parse_and_clean_generative_model( model, param );
    
    disp(upper('Adding documentation to model structure'));
    model = add_documentation_to_model( model, param );
    
    disp(upper('Adding parameters to model structure'));
    model = parse_and_clean_generative_model( model, param );
    
    disp(' '); print_simple_title('Clean up workspace and environment');
    model = clean_up_and_wrap_up_model( model, param );
    answer = clean_up( param );
end%3D
diary off
end%img2slml

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             HELPER METHODS                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% HELPER METHOD - CLEAN UP %%%%%%%%%%%%%%%%%%%%%%%%
function answer = clean_up( param )
%removes temp folder if debug set to true
if ~param.debug
    disp( 'Removing temporary folder' );
    if exist( [ pwd filesep 'temp' ] )
        rmdir( [ pwd filesep 'temp' ], 's'  );
    end
end

disp('Checking if model file exists on disk')
if exist( [param.model.filename(1:end-3) 'mat'] )
    answer = true;
else
    answer = false;
end

if isfield( param.protein, 'class' ) && ...
        strcmpi( param.protein.class, 'standardized_voxels' ) && ...
        exist( [param.model.filename(1:end-3) 'mat'] )
    delete( [param.model.filename(1:end-3) 'mat'] )
end
end%clean_up

function [dnaImagesDirectoryPath, cellImagesDirectoryPath, ...
    proteinImagesDirectoryPath, labels] = check_images_directory_paths( ...
    dnaImagesDirectoryPath, cellImagesDirectoryPath, ...
    proteinImagesDirectoryPath, param )

if ismember( param.train.flag, {'nuclear'} )
    if strcmpi( param.cell.class, 'framework' ) && ...
            strcmpi( param.cell.type, 'pca' ) && ...
            strcmpi( param.cell.class, 'framework' ) && ...
            strcmpi( param.cell.type, 'pca' )
        proteinImagesDirectoryPath = {};
    else
        cellImagesDirectoryPath = {};
        proteinImagesDirectoryPath = {};
    end
elseif ismember( param.train.flag, {'cell'} )
    if strcmpi( param.cell.class, 'framework' ) && ...
            strcmpi( param.cell.type, 'pca' ) && ...
            strcmpi( param.cell.class, 'framework' ) && ...
            strcmpi( param.cell.type, 'pca' )
        proteinImagesDirectoryPath = {};
    end
    
    if strcmpi( param.cell.class, 'cell_membrane' ) && ...
            strcmpi( param.cell.type, 'spharm_rpdm' )
        proteinImagesDirectoryPath = {};
        dnaImagesDirectoryPath = cellImagesDirectoryPath;
    end
    
    if strcmpi( param.cell.class, 'cell_membrane' ) && ...
            strcmpi( param.cell.type, 'ratio' )
        warning('Unable to train a cell membrane ratio model without a nuclear membrane');
        dnaImagesDirectoryPath = {};
        cellImagesDirectoryPath = {};
        proteinImagesDirectoryPath = {};
    end
elseif ismember( param.train.flag, {'framework'} )
    proteinImagesDirectoryPath = {};
elseif ismember( param.train.flag, {'protein'} )
    if isempty( proteinImagesDirectoryPath )
        warning('Unable to train a protein shape model without images');
        dnaImagesDirectoryPath = {};
        cellImagesDirectoryPath = {};
        proteinImagesDirectoryPath = {};
    end
else %param.train.flag == 'all'
    if isempty( proteinImagesDirectoryPath )
        warning('Unable to train a protein shape model without images');
        dnaImagesDirectoryPath = {};
        cellImagesDirectoryPath = {};
        proteinImagesDirectoryPath = {};
        labels = {};
    end
end

disp('Checking if using multiple datasets');
check_images_directory = @(x)( (isempty(x)) || (~isempty(x) && all(cellfun(@iscell,x))) );
if iscell( dnaImagesDirectoryPath ) && ...
        ( check_images_directory(dnaImagesDirectoryPath) && ...
        check_images_directory(cellImagesDirectoryPath) && ...
        check_images_directory(proteinImagesDirectoryPath) )
    disp('Multiple datasets found');
    disp('Checking consistency accross datasets')
    if numel(unique(nonzeros([length(dnaImagesDirectoryPath), ...
            length(cellImagesDirectoryPath), ...
            length(proteinImagesDirectoryPath)]))) == 1
        disp('All nonempty datasets have the same length')
    end
    
    labels = {};
elseif ~isempty(dnaImagesDirectoryPath)
    disp('Only one dataset found');
    labels = {};
else
    dnaImagesDirectoryPath = {};
    cellImagesDirectoryPath = {};
    proteinImagesDirectoryPath = {};
    labels = {};
end
end%check_images_directory_paths

function imgDir = parse_ometiff_deployed(arr)
file_array = {};
ch_num = {};
time_num = {};

disp('Parsing cell array')
for i = 1:(length(arr))
    split      = strsplit(arr{i},':');
    file_array = [file_array ml_ls(split{1})];
    
    %Check if given array has delimiter
    if length(strsplit(arr{1},':')) > 1
        ch_num     = [ch_num split{2}];
        if length(strsplit(arr{1},':')) > 2
            time_num   = [time_num split{3}];
        end
    end
end

%If no delimiter, no need to get ometiff func. handles
if isempty(ch_num)
    imgDir = file_array;
    return
end

imgDir = {};
channel_num = str2num(ch_num{1});
for i = 1:length(file_array)
    disp(['Parsing image ' file_array{i} ]);
    temp = get_list_of_function_handles_from_ometiff( [file_array{i}], channel_num);
    for j=1:1:length(temp)
        imgDir{length(imgDir)+1} = temp{j};
    end
end
end%parse_ometiff_deployed(arr)

%%%%%%%%%%%%%%% HELPER METHOD - CLEAN UP AND WRAP UP MODEL %%%%%%%%%%%%%%%%%
function model = clean_up_and_wrap_up_model( model, param )

filename = [param.model.filename(1:end-3) 'mat'];

if isdir( filename )
    if strcmpi( filename(end), filesep )
        filename = [ filename 'model.mat' ];
    else
        filename = [ filename filesep 'model.mat' ];
    end
else
    [ path, name, extension ] = fileparts( filename );
    
    if strcmpi( extension, '.mat' )
        filename = filename;
    else
        if isempty( path )
            filename = [ pwd filesep name '.mat' ];
        else
            filename = [ path filesep name '.mat' ];
        end
    end
end

if isfield( param, 'protein' ) && ...
        isfield( param.protein, 'class' ) && ...
        strcmpi( param.protein.class, 'standardized_voxels' ) && ...
        strcmpi( param.protein.type, 'standardized_map_half-ellipsoid' )
    if exist([pwd filesep 'model.mat'])
        delete([pwd filesep 'model.mat']);
    end
end

if ismember( param.train.flag, {'cell'} )
    if isfield( model, 'nuclearShapeModel' )
        model = rmfield( model, 'nuclearShapeModel' );
    end
    
    if isfield( model, 'proteinModel' )
        model = rmfield( model, 'proteinModel' );
    end
end

if ismember( param.train.flag, {'protein'} )
    if isfield( model, 'nuclearShapeModel' )
        model = rmfield( model, 'nuclearShapeModel' );
    end
    
    if isfield( model, 'cellShapeModel' )
        model = rmfield( model, 'cellShapeModel' );
    end
end

if ismember( param.train.flag, {'nuclear'} )
    if isfield( model, 'cellShapeModel' )
        model = rmfield( model, 'cellShapeModel' );
    end
    
    if isfield( model, 'proteinShapeModel' )
        model = rmfield( model, 'proteinShapeModel' );
    end
end

disp('Saving model structure to disk');
try
    save( filename, 'model', '-v7.3' );
catch
    save( filename, 'model' );
end
end%clean_up_and_wrap_up_model
