function model = parse_and_clean_generative_model( model, param )

try
    if ~isempty( param.model.name )
        model.name = param.model.name;
    end
catch
    if param.debug
        disp('Model name not set.')
    end
    model.name = 'unset';
end
disp(['model.name set to ''' model.name '''']);

try
    if ~isempty( param.model.id )
        model.id = param.model.id;
    end
catch
    if param.debug
        disp('Model ID not set.')
    end
    model.id = 'unset';
end
disp(['model.id set to ''' model.id'''']);

try
    if ~isempty( param.model.filename )
        model.filename = param.model.filename;
    end
catch
    if param.debug
        disp('Model filename not set.')
    end
    model.filename = 'unset';
end
disp(['model.filename set to ''' model.filename '''']);

if ismember( param.train.flag, {'nuclear', 'framework', 'all'} )
    disp(' '); print_comment_header('Adding parameters to model.nuclearShapeModel')
    try
        if ~isempty( param.nucleus.name )
            model.nuclearShapeModel.name = param.nucleus.name;
        end
    catch
        if param.debug
            disp('Nuclear shape model name not set.')
        end
        model.nuclearShapeModel.name = 'unset';
    end
    disp(['model.nuclearShapeModel.name set to ''' model.nuclearShapeModel.name '''']);

    try
        if ~isempty( param.nucleus.id )
            model.nuclearShapeModel.id = param.nucleus.id;
        end
    catch
        if param.debug
            disp('Nuclear shape model ID not set.')
        end
        model.nuclearShapeModel.id = 'unset';
    end
    disp(['model.nuclearShapeModel.id set to ''' model.nuclearShapeModel.id '''']);

    model.nuclearShapeModel.class = param.nucleus.class;
    disp(['model.nuclearShapeModel.class set to ''' model.nuclearShapeModel.class '''']);

    model.nuclearShapeModel.type = param.nucleus.type;
    disp(['model.nuclearShapeModel.type set to ''' model.nuclearShapeModel.type '''']);

    model.nuclearShapeModel.resolution = param.model.resolution;
end

disp(' '); print_comment_header('Adding parameters to model.cellShapeModel')
if ismember( param.train.flag, {'cell', 'framework', 'all'} )
    try
        if ~isempty( param.cell.name )
            model.cellShapeModel.name = param.cell.name;
        end
    catch
        if param.debug
            disp('Cell shape model name not set.')
        end
        model.cellShapeModel.name = 'unset';
    end
    disp(['model.cellShapeModel.name set to ''' model.cellShapeModel.name '''']);

    try
        if ~isempty( param.cell.id )
            model.cellShapeModel.id = param.cell.id;
        end
    catch
        if param.debug
            disp('Cell shape model ID not set.')
        end
        model.cellShapeModel.id = 'unset';
    end

    if ismember( param.train.flag, {'cell', 'framework', 'all' } )
	disp(['model.cellShapeModel.id set to ''' model.cellShapeModel.id '''']);
    	model.cellShapeModel.class = param.cell.class;
	disp(['model.nuclearShapeModel.class set to ''' model.cellShapeModel.class '''']);

    	model.cellShapeModel.type = param.cell.type;
    	disp(['model.cellShapeModel.type set to ''' model.cellShapeModel.type '''']);

    	model.cellShapeModel.resolution = param.model.resolution;
     end
end

disp(' '); print_comment_header('Adding parameters to model.proteinModel')
if ismember( param.train.flag, {'all', 'protein'} )
    try
        if ~isempty( param.protein.name )
            model.proteinModel.name = param.protein.name;
        end
    catch
        if param.debug
            disp('Protein shape model name not set.')
        end
        model.protein.name = 'unset';
    end
    disp(['model.proteinModel.name set to ''' model.proteinModel.name '''']);

    % protein.id  (optional) Holds the id of the protein model. The default is empty.
    try
        if ~isempty( param.protein.id )
            model.proteinModel.id = param.protein.id;
        end
    catch
        if param.debug
            disp('Protein shape model ID not set.')
        end
        model.proteinModel.id = 'unset';
    end
    disp(['model.proteinModel.id set to ''' model.proteinModel.id '''']);

    model.proteinModel.class = param.protein.class;
    disp(['model.proteinModel.class set to ''' model.proteinModel.class '''']);

    model.proteinModel.type = param.protein.type;
    disp(['model.proteinModel.type set to ''' model.proteinModel.type '''']);

    model.proteinModel.resolution = param.model.resolution;
else
    disp('Model is empty. Nothing to do.');
end

parameterization = [ pwd filesep 'param' ];
files = dir( [parameterization filesep 'param*.*'] );
model.dataset.parameterization = {};
if isfield(param,'save_segmentations') && param.save_segmentations == true 
    model.dataset.segmentation = {};
end

for i=1:1:length(files)
    model.dataset.parameterization{length(model.dataset.parameterization)+1} = ...
        files(i).name;
    if isfield(model.dataset,'segmentation') && exist([parameterization filesep 'param' num2str(i) filesep 'seg.mat'])
        seg = load([parameterization filesep 'param' num2str(i) filesep 'seg.mat']);
        seg = seg.seg;
        if isfield(seg, 'nuc') && isfield(seg, 'cell')
            model.dataset.segmentation{length(model.dataset.segmentation)+1} = ...
                struct('nucleus',seg.nuc,'cell',seg.cell);
        elseif ~isfield(seg, 'nuc') && isfield(seg, 'cell')
            model.dataset.segmentation{length(model.dataset.segmentation)+1} = ...
                struct('cell',seg.cell);
        else
            model.dataset.segmentation{length(model.dataset.segmentation)+1} = ...
                struct('nucleus',seg.nuc);
        end
            
    elseif isfield(model.dataset, 'segmentation')
        model.dataset.segmentation{length(model.dataset.segmentation)+1} = 'Empty';
    end
end
end%parse_and_clean_generative_model
