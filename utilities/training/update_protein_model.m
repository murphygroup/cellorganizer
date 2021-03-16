%%%%%%%%%%%%%%%%%%% HELPER METHOD - UPDATE PROTEIN MODEL %%%%%%%%%%%%%%%%%%
function model = update_protein_model( dimensionality, model, param )

if strcmpi( dimensionality, '3D' )
    % protein.id  (optional) Holds the id of the protein model. The default is empty.
    try
        model.proteinModel.id = param.protein.id;
    catch
        if param.debug
            disp('Protein shape model ID not set. Ignoring meta data field.')
        end
    end
    
    % protein.type              (optional) Holds the protein model type. The default is "vesicle".
    try
        model.proteinModel.type = param.protein.type;
    catch
        if param.debug
            disp('Protein shape model type not set. Using default: vesicle.')
            model.proteinModelModel.type = 'vesicle';
        end
    end
    
    % protein.class             Holds the protein class, e.g. lysosome, endosome.
    try
        model.proteinModel.class = param.protein.class;
    catch
        if param.debug
            disp('Protein shape model class not set. Ignoring meta data field.')
        end
    end
    
    % protein.cytonuclearflag   (optional) Determines whether the protein pattern will be generated in
    %                           the cytosolic space ('cyto'), nuclear space ('nuc') or everywhere ('all').
    %                           Default is cyto.
    try
        flag = param.protein.cytonuclearflag;
        if ~strcmpi( flag, 'cyto' ) && ~strcmpi( flag, 'nuc' ) ...
                ... & ~strcmpi( flag, 'all' )
                if param.debug
                warning('Protein shape model cytonuclear flag is unrecognized. Using default: cyto.');
                end
                flag = 'cyto';
        end
        
        model.proteinModel.cytonuclearflag = flag;
        clear flag;
    catch
        if param.debug
            disp( 'Protein shape model cytonuclear flag not set. Using default: cyto.');
        end
        model.proteinModel.cytonuclearflag = 'cyto';
    end
    
    %set model dimensionality
    model.proteinModel.dimensionality = dimensionality;
    
    %dpsulliv 2/24/13
    %There are no assumptions of similar resolution, therefore there
    %must be separate resolutions tracked for framework and protein
    model.nuclearShapeModel.resolution = param.model.resolution;
    model.cellShapeModel.resolution = param.model.resolution;
    
    %icaoberg 4/23/2014
    %included missing protein shape resolution;
    %fixed bug where it was saving the wrong resolution. it was saving
    %the downsampled resolution instead of the protein resolution
    model.proteinModel.resolution = param.model.protein_resolution;
end
end%update_protein_model
