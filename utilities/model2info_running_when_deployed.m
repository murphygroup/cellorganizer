function model2info_running_when_deployed( model, fileID, options )
% MODEL2INFO_RUNNING_WHEN_DEPLOYED List basic information about an SLML model.
%
% model   model structure (required)
% options (optional) options for modifying information displayed

% Ivan E. Cao-Berg
%
% Copyright (C) 2020 Murphy Lab
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
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

% 2018/03/18 icaoberg Updated method to save results to a file

%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL IDENTIFIERS %%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID, '## Model Identifiers\n');
fprintf( fileID,'```\n' );
identifiers = {'filename', 'name', 'id'};
for j = 1:3
    identifier = identifiers{j};
    if isfield(model, identifier)
        identifier_value = eval(sprintf('model.%s', identifier));
        fprintf(fileID, 'model.%s = ''%s'';\n', identifier, identifier_value);
    else
        fprintf(fileID, 'model.%s = ''%s'';\n', identifier, 'UNSET');
    end
end
fprintf( fileID,'```\n' );

%%%%%%%%%%%%%%%%%%%%%%%%%% DIMENSIONALITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( fileID,'## Dimensionality\n' );
fprintf( fileID,'```\n' );
if isfield( model, 'dimensionality' )
    fprintf( fileID, 'model.dimensionality = ''%s'';\n', model.dimensionality );
else
    fprintf(fileID, 'model.dimensionality = ''%s'';\n', model.dimensionality, 'UNSET');
    warning( 'This is an invalid model' );
end
fprintf( fileID,'```\n' );

%%%%%%%%%%%%%%%%%%%%%%%%%% DOCUMENTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'documentation' )
    fprintf( fileID,'## Documentation\n' );
    fprintf( fileID,'```\n' );
    fnames = fieldnames( model.documentation );
    for i=1:1:length(fnames)
        fname = fnames{i};
        fvalue = eval(sprintf('model.documentation.%s', fname));
        fprintf(fileID, 'model.documentation.%s = ''%s'';\n', fname, fvalue);
        clear fname
    end
    fprintf( fileID,'```\n' );
    clear fnames
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBMODEL INFORMATION %%%%%%%%%%%%%%%%%%%%%%%
fprintf( fileID,'## Nuclear shape model information\n' );
fprintf( fileID,'```\n' );
if isfield(model, 'nuclearShapeModel')
    model2info2(model.nuclearShapeModel, 'nuclearShapeModel', fileID)
else
    fprintf( fileID,'This model does not contain a nuclear shape submodel.\n' );
end
fprintf( fileID,'```\n' );

fprintf( fileID,'## Cell shape model information\n' );
fprintf( fileID,'```\n' );
if isfield(model, 'cellShapeModel')
    model2info2(model.cellShapeModel, 'cellShapeModel', fileID)
else
    fprintf( fileID,'This model does not contain a cell shape submodel.\n' );
end
fprintf( fileID,'```\n' );

fprintf( fileID,'## Protein pattern model information\n' );
fprintf( fileID,'```\n' );
if isfield(model, 'proteinModel')
    model2info2(model.proteinModel, 'proteinModel', fileID)
else
    fprintf( fileID,'This model does not contain a protein pattern submodel.\n' );
end
fprintf( fileID,'```\n' );

fprintf( fileID,'\n' );
disp('Checking Type of Model')
%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATASET Information %%%%%%%%%%%%%%%%%%%%%%%%%%
% if isfield(model, 'dataset')
% 	answer = datasetinfo2html(model, fileID);
% 	if isfield(model.dataset, 'segmentation')
% 		ans_seg = segmentation2html(model, fileID, 'report');
%     end
% end

%   if isfield(options, 'labels')
%       labels = options.labels;
%   else
%       labels = [];
%   end


%%%%%%%%%%%%%%%%%%%%%%%% IS PPM? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_ppm_model( model )
    fprintf( fileID,'## PPM model information \n' );
    fprintf( fileID,'```\n' );
    if isfield(model, 'proteinModel.ppm')
        model2info3 (model.proteinModel.ppm, 'proteinModel.ppm', fileID)
    else
        fprintf( fileID,'This model does not contain a PPM submodel.\n' );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%% IS DIFFEOMORPHIC? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_diffeomorphic( model )
    fprintf(fileID,'## Additional Figures\n');
    f = figure('visible','off');
    nimgs = size(model.cellShapeModel.positions,1);
    nlabels = size(model.cellShapeModel.positions,1);
    labels = reshape(repmat([1:nlabels],ceil(nimgs/nlabels),1),[],1);
    options.plot_dims = [1,2];
    options.subsize = 1000;
    showShapeSpaceFigure( model, labels, options );
    saveas( f, 'show_shape_space.png', 'png' );
    I = imread( 'show_shape_space.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'show_shape_space_thumbnail.png' );
    fprintf( fileID, '[![Shape space](./show_shape_space_thumbnail.png)](./show_shape_space.png)\n' );
    movefile( '*.png', './report/' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS SPHARM-obj? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_spharm_obj_model( model )
    fprintf( fileID,'## SPHARM-obj model information \n' );
    fprintf( fileID,'```\n' );
    if isfield(model, 'proteinModel.spharm_obj_model.cellShapeModel')
        model2info2(model.proteinModel.spharm_obj_model.cellShapeModel, ...
            'proteinModel.spharm_obj_model.cellShapeModel', fileID)
    else
        fprintf( fileID,'This model does not contain a SPHARM-obj submodel.\n' );
    end
    fprintf( fileID,'```\n' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS TCell MODEL? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_tcell_model( model )
    fprintf(fileID,'%%%%%% Additional Figures\n');
    ShowTcellEnrichment(model.model_filename, options);
    img_name=ml_ls('*model_mean_enrichment_plot_mdl*.png');
    I = imread( img_name{1} );
    imwrite( I, 'show_enrichment.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'show_enrichment_thumbnail.png' );
    
    fprintf( fileID, '%%\n' );
    fprintf( fileID, '%% <html>\n' );
    fprintf( fileID, '%% <figure>\n' );
    fprintf( fileID, '%% <a href="show_enrichment.png" ><img src="show_enrichment_thumbnail.png" /></a>\n' );
    fprintf( fileID, '%% <figcaption>Figure 1. Shape enrichment for input model. </figcaption>\n' );
    fprintf( fileID, '%% </figure>\n' );
    fprintf( fileID, '%% </html>\n' );
    fprintf( fileID, '%%\n' );
    close all
    
    movefile( '*.png', './report' );
end

%%%%%%%%%%%%%%%%%%%%%%%% IS PCA FRAMEWORK? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_pca_framework( model )
    fprintf(fileID,'\n## Additional Figures\n');
    f = figure('visible','off');
    if isfield(options, 'labels')
        labels = options.labels;
    else
        labels = [];
    end
    showPCAShapeSpaceFigure(model, labels, options);
    saveas( f, 'show_shape_space.png', 'png' );
    I = imread( 'show_shape_space.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'show_shape_space_thumbnail.png' );
    fprintf( fileID, '[![Shape space](./show_shape_space_thumbnail.png)](./show_shape_space.png)\n' );
    movefile( '*.png', './report/' );
end

%%%%%%%%%%%%%%%%%%%%%%%% IS CLASSIC MODEL? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_classic_model( model )
    fprintf(fileID,'\n## Additional Figures\n');
    f = figure('visible','off');
    fprintf( fileID,'There are no additional figures for this model at the moment.\n' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS SPHARM MODEL? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isfield(model,'proteinModel.spharm_obj_model'))
    model=model.proteinModel.spharm_obj_model;
end

if is_spharm_model( model )
    if(isfield(model,'proteinModel'))
        model=model.proteinModel.spharm_obj_model;
    end
    disp('Writing figures for SPHARM MODEL')
    fprintf(fileID,'%%%%%% Additional Figures\n');
    f = figure('visible','off');
    
    if(strcmp(options.shape_evolution,'null')==0)
        spharm_rpdm_sample_or_reconstruct_images_figure(model,options);
        saveas( gcf, 'show_shape_evolution.png', 'png' );
        I = imread( 'show_shape_evolution.png' );
        I = imresize( I, 0.50 );
        imwrite( I, 'show_shape_evolution_thumbnail.png' );
        fprintf( fileID, '%% <figure>\n' );
        fprintf( fileID, '%% <a href="show_shape_evolution.png"><img src="show_shape_evolution_thumbnail.png"></a>\n' );
        fprintf( fileID, '%% <figcaption>Shape evolution. </figcaption>\n' );
        fprintf( fileID, '%% </figure>\n' );
        close all
        f = figure('visible','off');
    end
    %     show_SPHARM_RPDM_shape_evolution_figure(model);
    %     saveas( gcf, 'show_shape_evolution.png', 'png' );
    %     I = imread( 'show_shape_evolution.png' );
    %     I = imresize( I, 0.50 );
    %     imwrite( I, 'show_shape_evolution_thumbnail.png' );
    %     fprintf( fileID, '%%\n' );
    %     fprintf( fileID, '%% <html>\n' );
    %     fprintf( fileID, '%% <figure>\n' );
    %     fprintf( fileID, '%% <a href="show_shape_evolution.png" ><img src="show_shape_evolution_thumbnail.png" /></a>\n' );
    %     fprintf( fileID, '%% <figcaption>Figure 1. Shape evolution. </figcaption>\n' );
    %     fprintf( fileID, '%% </figure>\n' );
    %     close all
    %     f = figure('visible','off');
    show_SPHARM_RPDM_Shape_Space_Figure(model);
    saveas( f, 'show_shape_space.png', 'png' );
    I = imread( 'show_shape_space.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'show_shape_space_thumbnail.png' );
    fprintf( fileID, '%% <figure>\n' );
    fprintf( fileID, '%% <a href="show_shape_space.png" ><img src="show_shape_space_thumbnail.png" /></a>\n' );
    fprintf( fileID, '%% <figcaption>Figure 2. Shape space. </figcaption>\n' );
    fprintf( fileID, '%% </figure>\n' );
    fprintf( fileID, '%% </html>\n' );
    fprintf( fileID, '%%\n' );
    close all
    
    movefile( '*.png', './report' );
end
end

function model2info2( model, submodel, fileID )
attributes = {'name', 'id', 'class', 'type', 'resolution'};
for j = 1:length(attributes)
    attribute = attributes{j};
    if isfield(model, attribute)
        attribute_value = eval(sprintf('model.%s', attribute));
        switch class(attribute_value)
            case 'double'
                if length(attribute_value) == 1
                    fprintf(fileID, 'model.%s.%s = %d;\n', submodel, attribute, attribute_value);
                else
                    attribute_value_string = ['[' num2str(attribute_value(1))];
                    for i=2:length(attribute_value)
                        attribute_value_string = [attribute_value_string ' ' num2str(attribute_value(i))];
                    end
                    attribute_value_string = [attribute_value_string ']'];
                    fprintf(fileID, 'model.%s.%s = %s;\n', submodel, attribute, attribute_value_string);
                end
            otherwise
                fprintf(fileID, 'model.%s.%s = ''%s'';\n', submodel, attribute, attribute_value);
        end
    else
        fprintf(fileID, 'model.%s.%s = ''%s'';\n', submodel, attribute, 'UNSET');
    end
end
end%model2info2

function model2info3( model, submodel, fileID )
attributes = {'cellmembrane', 'nuclear', 'full'};
for k = 1:length(attributes)
    attribute = attributes{k};
    if isfield(model, attribute)
        attribute_value = eval(sprintf('model.%s', attribute));
        switch class(attribute_value)
            case 'double'
                if length(attribute_value) == 1
                    fprintf(fileID, 'model.%s.%s = %d;\n', submodel, attribute, attribute_value);
                else
                    attribute_value_string = ['[' num2str(attribute_value(1))];
                    for i=2:length(attribute_value)
                        attribute_value_string = [attribute_value_string ' ' num2str(attribute_value(i))];
                    end
                    attribute_value_string = [attribute_value_string ']'];
                    fprintf(fileID, 'model.%s.%s = %s;\n', submodel, attribute, attribute_value_string);
                end
            otherwise
                fprintf(fileID, 'model.%s.%s = ''%s'';\n', submodel, attribute, attribute_value);
        end
    else
        fprintf(fileID, 'model.%s.%s = ''%s'';\n', submodel, attribute, 'UNSET');
    end
end
end%model2info3
