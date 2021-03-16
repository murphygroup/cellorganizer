%%%%%%%%%%%%%%%%% HELPER METHOD - CLEAN UP INPUT ARGUMENTS %%%%%%%%%%%%%%%%
function param = clean_up_training_input_arguments( dimensionality, param )

disp(' '); print_simple_title('Checking model dimensionality' );
if ~(strcmpi( dimensionality, '2d' ) || ...
        strcmpi( dimensionality, '3d' ))
    warning('Unrecognized or unsupported dimensionality. Exiting method.');
    param = [];
    return
else
    %icaoberg 9/16/2013
    param.model.dimensionality = dimensionality;
end
disp(['Dimensionality set to ''' dimensionality '''']);

disp(' '); print_simple_title('Checking the contents of options structure' );
disp([upper('Checking') ' param.model.name'])
try
    temp = param.model.name;
    if ~isa( temp, 'char' )
        warning(['Optional input argument param.model.name must be a string. Changing to default value']);
        param.model.name = 'unset';
    end 
catch
    param.model.name = 'unset';
end
disp(['param.model.name set to ''' param.model.name '''']);

disp([upper('Checking') ' param.model.id'])
try
    temp = param.model.id;
    %icaoberg 7/9/2012
    if ~isa( temp, 'char' )
        warning(['CellOrganizer: Optional input argument param.model.id must be a string. Changing to default value']);
        param.model.id = uuidgen();
    end
catch
    param.model.id = uuidgen();
end
disp(['param.model.id set to ''' param.model.id '''']);

disp([upper('Checking') ' param.downsampling'])
try
    downsampling = param.downsampling;
catch
    if strcmpi( dimensionality, '3D' )
        downsampling = [5 5 1];
    else
        downsampling = [1 1];
    end
end
param.downsampling = downsampling;
disp(['param.downsampling set to [' num2str(param.downsampling) ']']);

disp([upper('Checking') ' param.model.filename'])
try
    temp = param.model.filename;
    if ~isa( temp, 'char' )
        warning('Input parameter filename must be a string.');
        param.model.filename = 'model.xml';
    end
catch err
    warning( 'Unable to set output filename. Setting default value.' );
    param.model.filename = 'model.xml';
end
disp(['param.model.filename set to ''' param.model.filename '''']);

disp('Checking param.model.resolution')
try
    disp('Setting param.model.original_resolution to param.model.resolution');
    param.model.original_resolution = param.model.resolution;
    %D. Sullivan 2/24/13 added param.model.protein_resolution
    disp('Setting param.model.protein_resolution to param.model.resolution');
    param.model.protein_resolution = param.model.original_resolution;
    disp('Setting downsampling vector');
    param.downsampling = downsampling;
    param.model.resolution_orig = param.model.resolution;
    disp('Computing param.model.resolution')
    param.model.resolution = param.model.resolution .* downsampling;
    clear downsampling
catch
    disp('Model resolution not set. Exiting method.');
    param = [];
    return
end

disp(' '); print_simple_title('Checking selection of model class/type' );
%%%%%%%%%% CHECKING SELECTION OF NUCLEAR MEMBRANE MODEL CLASS/TYPE %%%%
disp(upper('CHECKING SELECTION OF NUCLEAR MEMBRANE MODEL CLASS/TYPE'))
if ismember( param.train.flag, {'nuclear','framework', 'all'} )
    disp('Checking selection of nuclear model class/type' );
    answer = selection_flag( dimensionality, 'nucleus', param.nucleus.class, param.nucleus.type );
    if ~answer
        param = [];
        return
    else
        disp(['Nuclear shape model class/type:' param.nucleus.class '/' param.nucleus.type ' is valid']);
    end
end

%%%%%%%%%%% CHECKING SELECTION OF CELL MEMBRANE MODEL CLASS/TYPE %%%%%%
disp(upper('CHECKING SELECTION OF CELL MEMBRANE MODEL CLASS/TYPE'))
if ismember( param.train.flag, {'cell','framework', 'all'} )
    disp('Checking selection of cell membrane model class/type' );
    answer = selection_flag( dimensionality, 'cell', param.cell.class, param.cell.type );
    if ~answer
        param = [];
        return
    else
        disp(['Cell membrane shape model class/type:' param.cell.class '/' param.cell.type ' is valid']);
    end
end

%%%%%%%%%% CHECKING SELECTION OF NUCLEAR MEMBRANE MODEL CLASS/TYPE %%%%
disp(upper('CHECKING SELECTION OF NUCLEAR MEMBRANE MODEL CLASS/TYPE'))
if ismember( param.train.flag, {'protein','all'} )
    disp('Checking selection of protein model class/type' );
    answer = selection_flag( dimensionality, 'protein', param.protein.class, param.protein.type );
    if ~answer
        param = [];
        return
    else
        disp(['Protein model class/type:' param.protein.class '/' param.protein.type ' is valid']);
    end
end

disp(' '); print_simple_title('Checking contents of other options' );
disp([upper('Checking the contents of ') ...
    'param.nucleus, param.cell ' upper('and') ' param.protein']);
compartments = { 'nucleus', 'cell', 'protein' };
fields = { 'type', 'name', 'id' };
for c=1:1:length(compartments)
    compartment = [ 'param.' compartments{c} ];
    for f=1:1:length(fields)
        field = [ compartment '.' fields{f} ];
        try
            eval(['temp=' field ';']);
            if ~isa( temp, 'char' )
                eval([field '=''''' ';']);
            end
        catch
            eval([field '=''''' ';']);
        end
    end
end

disp([upper('Checking the contents of ') 'param.documentation']);
try
    documentation = param.documentation;
    if ~isa( documentation, 'struct' )
        error('Documentation must be a structure.');
    end
catch
    param.documentation = [];
end
if isempty( param.documentation )
    disp('param.documentation if empty');
else
    disp('param.documentation set and not empty.')
end

disp([upper('Checking the contents of ') 'param.train.flag']);
try
    trainFlag = param.train.flag;
    if ~isa( trainFlag, 'char' )
        error('Training flag must be a string');
    end
catch
    param.train.flag = 'all';
    %icaoberg july 5, 2012
    trainFlag = param.train.flag;
end
disp(['param.train.flag set to ''' param.train.flag '''']);

end%clean_up_training_input_arguments
