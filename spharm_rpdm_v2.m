function answer = spharm_rpdm_v2(spharm_obj_files,options)
%call to spharm_rpdm
% 1/26/2021 R.F. Murphy - save indices of objects that were successfully parameterized
% 2/8/2021 R.F. Murphy - set a default hd_thresh
% 1/31/2023 R.F. Murphy - correct passing of objects (was using nuc_path; s/b cell_path)
%                       add check for correct train.flag
% 4/25/2023 R.F.Murphy - change default for hd_thresh to spharm_rpdm.hd_thresh
%                       for compatibility with other functions

%IMG2SLML2
answer = false;

if ~exist( [ pwd filesep 'log'], 'dir' )
    mkdir( 'log' );
end
c=clock;
logfile = '';

for i=1:1:length(c)
    logfile = ['',logfile,num2str(c(i))]; %#ok<AGROW>
end

logfile = [ pwd filesep 'log' filesep logfile, '.log' ];
diary( logfile )

disp('Running img2slml')
disp('Check number of input arguments')
dnaImagesDirectoryPath = spharm_obj_files;
cellImagesDirectoryPath = spharm_obj_files;

 disp('Checking and getting default parameters')
 options = get_cellorganizer_default_parameters( 'training', options );
 
 [spharm_obj_files, spharm_obj_files, ...
  ~,labels] = check_images_directory_paths( ...
        dnaImagesDirectoryPath, cellImagesDirectoryPath, ...
        [], options );
 options = clean_up_training_input_arguments( '3D', options );
 
 if options.verbose
    fprintf( 1, '%s\n', 'Checking the existence of temporary folder' );
end
 
 % %TRAIN THE GENERATIVE MODEL
disp(' '); print_large_title('Training generative model' );
options.dimensionality = '3D';

%IMG2SLML2


%IMG2MODEL2

model = [];

default_options = struct( 'display', false ,...
    'debug', false, ...
    'verbose', false, ...
    'paramdir', [pwd filesep 'param'], ...
    'tempparent', [pwd filesep 'temp'])
default_options.spharm_rpdm = struct('hd_thresh', 20);
options = ml_initparam( options, default_options);

check_required_options(options, {'dimensionality'})
check_required_options(options.model, {'resolution'})

disp('Setting up data');

%Could be cleaned up
[spharm_obj_files, spharm_obj_files,~, options] = setup_data(spharm_obj_files, ...
    spharm_obj_files, {}, options);
%Could be cleaned up 

%Could be cleaned up
disp('Setting up model options');
options = setup_model_options(spharm_obj_files, spharm_obj_files, {}, options);
%Could be cleaned up

% check if train flag option matches requirements
if options.train.flag ~= 'cell'
    warning(['options.train.flag was set to ' options.train.flag '; should be "cell"']);
    options.train.flag = 'cell';
end

paramfiles = cell(size(spharm_obj_files));
isdone = false(size(spharm_obj_files));

disp(' '); print_large_title('Processing images');
for i = 1:length(spharm_obj_files)

    paramfiles{i} = [options.paramdir filesep 'param' num2str(i) '.mat'];
    
    fname = paramfiles{i};
    [fname_path, fname_name, fname_ext] = fileparts(fname);
    fname = [fname_path, fname_name];
    [can_start, final_name, final_exists, tmpfile] = chunk_start(fname, '.mat');
    if final_exists
        isdone(i) = true;
    end
    if ~can_start
        continue
    end
    
    [spharm_obj,options.cell_image_path] = readfileifnonblank(spharm_obj_files,i);
    %options.nuc_image_path = options.cell_image_path;
    [immask,options.crop_image_path] = readfileifnonblank(options.masks,i);

    
    if ~isa(spharm_obj, 'uint8' )
        spharm_obj = uint8(spharm_obj);
    end
    if ~isa(immask, 'uint8' )
        immask = uint8(immask);
    end
    
    savedir = [options.paramdir filesep 'param' num2str(i)];

    try
        [cell_params] = img2param(spharm_obj, spharm_obj, spharm_obj, immask, savedir, options);
        %Serena 03/21 - eliminate NaN values
        if ~isempty(cell_params)
            if (sum(isnan(cell_params.cell.fvec(:)))>0 || sum(isnan(cell_params.cell.vertices(:)))>0 || sum(isnan(cell_params.cell.faces(:)))>0 || sum(isnan(cell_params.cell.sph_verts(:)))>0)
                isdone(i)=false;
            else
                save(paramfiles{i}, '-struct', 'cell_params')
                isdone(i) = true;
            end
            chunk_finish(fname);
        end
    catch the_error
        warning(['Unable to extract parameters for cell ' num2str(i) ...
            ': it will be ignored.']);
        getReport( the_error )
        disp( 'Check the images exist or that you are using the correct options.' );
    end
end
goodparamfiles = paramfiles(isdone);

model=param2model(goodparamfiles,options);


if isempty( model )
    warning( ['Method img2model returned an empty model. Exiting method.']);
    answer = false;
    return
end

%
disp('in spharm_rpdm_v2: save indices of good objects to the model file');
%
model.cellShapeModel.parameterization_successful = isdone;

%PARSE GENERATIVE MODEL INTO SLML INSTANCE
disp(' '); print_large_title('Parse and clean generative model' );

disp(upper('Adding parameters to model structure'));
model = parse_and_clean_generative_model( model, options );
    
disp(upper('Adding documentation to model structure'));
model = add_documentation_to_model( model, options );
    
disp(upper('Adding parameters to model structure'));
model = parse_and_clean_generative_model( model, options );
    
disp(' '); print_simple_title('Clean up workspace and environment');
model = clean_up_and_wrap_up_model( model, options );
answer = clean_up( options );

diary off


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             HELPER METHODS                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_required_options(options, required_fields)
ismissing = false(size(required_fields));
for j = 1:length(required_fields)
    if ~isfield(options, required_fields{j})
        ismissing(j) = true;
    end
end
if any(ismissing)
    error(['Missing options fields: ' strjoin(required_fields(ismissing), ', ')])
end
end

function [dna_images_list, cell_images_list, protein_images_list, options] = setup_data(dna_images, cell_images, prot_images, options)

options = ml_initparam(options, struct('masks', []));

if ~exist(options.paramdir, 'dir')
    mkdir(options.paramdir);
end

disp(' '); print_simple_title('Creating list of nuclear membrane images');
dna_images_list = {};
label = -1;
temp_labels = {};
if isa(dna_images, 'cell' ) && ...
        ~any(cellfun( @(x)(isa(x,'function_handle')), dna_images ))
    for i=1:1:length((dna_images))
        dataset = dna_images{i};
        temp = ml_ls( dataset );
        if length(temp) ~= 1
            label = label+1;
        end
        if isempty(dna_images_list)
            disp('Adding first dataset to list');
            dna_images_list = temp;
            for j=1:1:length(dna_images_list)
                disp(['Adding file ' dna_images_list{j}]);
                temp_labels{end+1} = num2str(label);
            end
        else
            disp('Adding another dataset to list')
            for j=1:1:length(temp)
                dna_images_list{end+1} = temp{j};
                disp(['Adding file ' temp{j}]);
                temp_labels{end+1} = num2str(label);
            end
        end
    end
    if isempty(options.labels)
        options.labels = temp_labels;
    end
    clear temp_labels;
else
    if iscell( dna_images ) && ... 
        all(cellfun( @(x)(isa(x,'function_handle')), dna_images ))
        dataset = dna_images;
        label = label + 1;
        disp('Adding datasets to list')
        for j=1:1:length(dataset)
            dna_images_list{end+1} = dataset{j};
            disp(['Adding file function handle ' num2str(j) ' to list']);
            temp_labels{end+1} = num2str(label);
        end
        options.labels = temp_labels;
        clear dataset
    else
        dataset = dna_images;
        temp = ml_ls( dataset );
        label = label + 1;
        disp('Adding datasets to list')
        for j=1:1:length(temp)
            dna_images_list{end+1} = temp{j};
            disp(['Adding file ' temp{j}]);
            temp_labels{end+1} = num2str(label);
        end
        options.labels = temp_labels;
        clear temp_labels;
    end
end

if isempty(dna_images_list)
    disp('List is empty or no nuclear membrane images found');
end

disp(' '); print_simple_title('Creating list of cell membrane images');
cell_images_list = {};
if isa(cell_images, 'cell' ) && ...
        ~any(cellfun( @(x)(isa(x,'function_handle')), cell_images ))
    for l=1:1:length((cell_images))
        dataset = cell_images{l};
        temp = ml_ls( dataset );
        
        if isempty(cell_images_list)
            disp('Adding first dataset to list');
            cell_images_list = temp;
            for j=1:1:length(cell_images_list)
                disp(['Adding file ' cell_images_list{j}]);
            end
        else
            disp('Adding another dataset to list')
            for j=1:1:length(temp)
                cell_images_list{end+1} = temp{j};
                disp(['Adding file ' temp{j}]);
            end
        end
    end
else
    if iscell( cell_images ) && ...
            all(cellfun( @(x)(isa(x,'function_handle')), cell_images ))
        dataset = cell_images;
        label = label + 1;
        disp('Adding datasets to list')
        for j=1:1:length(dataset)
            cell_images_list{end+1} = dataset{j};
            disp(['Adding file function handle ' num2str(j) ' to list']);
        end
        clear dataset
    else
        dataset = cell_images;
        temp = ml_ls( dataset );
        disp('Adding datasets to list')
        for j=1:1:length(temp)
            cell_images_list{end+1} = temp{j};
            disp(['Adding file ' temp{j}]);
        end
    end
end

if isempty(cell_images_list)
    disp('List is empty or no cell membrane images found');
    cell_images_list = cell(size(dna_images_list));
end

disp(' '); print_simple_title('Creating list of protein pattern images');
protein_images_list = {};
if isa(prot_images, 'cell' ) && ...
        ~any(cellfun( @(x)(isa(x,'function_handle')), prot_images ))
    for i=1:1:length((prot_images))
        dataset = prot_images{i};
        temp = ml_ls( dataset );
        
        if isempty(protein_images_list)
            disp('Adding first dataset to list');
            protein_images_list = temp;
            for j=1:1:length(protein_images_list)
                disp(['Adding file ' protein_images_list{j}]);
            end
        else
            disp('Adding another dataset to list')
            for j=1:1:length(temp)
                protein_images_list{end+1} = temp{j};
                disp(['Adding file ' temp{j}]);
            end
        end
    end
else
    if iscell( prot_images ) && ... 
            all(cellfun( @(x)(isa(x,'function_handle')), prot_images ))
        dataset = prot_images;
        label = label + 1;
        disp('Adding datasets to list')
        for j=1:1:length(dataset)
            protein_images_list{end+1} = dataset{j};
            disp(['Adding file function handle ' num2str(j) ' to list']);
            temp_labels{end+1} = num2str(label);
        end
        clear dataset
    else
        dataset = prot_images;
        temp = ml_ls( dataset );
        disp('Adding datasets to list')
        for j=1:1:length(temp)
            protein_images_list{end+1} = temp{j};
            disp(['Adding file ' temp{j}]);
        end
    end
end

if isempty(protein_images_list)
    disp('List is empty or no protein images found');
    protein_images_list = cell(size(dna_images_list));
end

if ischar(options.masks)
    options.masks = ml_ls(options.masks);
end
if isempty(options.masks)
    options.masks = cell(size(dna_images_list));
end

if isempty(cell_images_list) && ~isempty(options.protein.type) && ...
        strcmp(options.protein.type, 'standardized_map_half-ellipsoid')
    disp('Setting standardized_map_half-ellipsoid model' );
    [cell_images_list, options] = tcell_setup_options(options);
end

disp(' ' ); disp('Saving dataset and label information')
options.dataset.nuclear_membrane_images = dna_images_list;
options.dataset.cell_membrane_images = cell_images_list;
options.dataset.protein_images = protein_images_list;
options.dataset.labels = options.labels;
if isfield( options, 'labels' )
    options = rmfield( options, 'labels' );
end

model.dataset = options.dataset;


end

function options = setup_model_options(dna_images, cell_images, prot_images, options)
%this function sets the default model options for cellorganizer
component_struct = struct('type', '', ...
    'name', '', ...
    'id', '');

options = ml_initparam(options, struct('nucleus', []));
% xruan 01/05/2016 change ml_initparam(options, component_struct); to ml_initparam(options.nucleus, component_struct);
options.nucleus = ml_initparam(options.nucleus, component_struct);
options.documentation.numimgs = numel(dna_images);

if strcmpi(options.dimensionality, '2D')
    if isempty(options.nucleus.type)
        options.nucleus.type = 'medial axis';
    end
    
    if ~all(cellfun(@isempty, cell_images)) && isempty(options.cell.type)
        options.cell.type = 'ratio';
    end
    
    if ~all(cellfun(@isempty, prot_images)) && isempty(options.protein.type)
        options.protein.type = 'vesicle';
    end
    
elseif strcmpi(options.dimensionality, '3D')
    if isempty(options.nucleus.type)
        options.nucleus.type = 'cylindrical_surface';
    end
    
    if ~all(cellfun(@isempty, cell_images)) && isempty(options.cell.type)
        options.cell.type = 'ratio';
    end
    
    if ~all(cellfun(@isempty, prot_images)) && isempty(options.protein.type)
        options.protein.type = 'vesicle';
    end
else
    error('Unsupported dimensionality. Exiting method.')
end
end

function [img,filename] = readfileifnonblank(files,i)
if ~isempty(files)
    filename = files{i};
    
    if ~isempty( filename )
        if ~strcmpi(class(filename), 'function_handle')
            disp(['Reading file ' filename] );
        end
    end
    img = ml_readimage(filename);
else
    filename = [];
    img = [];
end
end

function filename = copypathifnonblank(files,i)
if ~isempty(files)
    filename = files{i};
else
    filename = [];
end
end

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
    disp('Checking consistency across datasets')
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
        if length(strsplit(arr{1},':')) > 1
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
        imgDir{end+1} = temp{j};
    end
end
end
