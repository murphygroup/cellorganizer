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
