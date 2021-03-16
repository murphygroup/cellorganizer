function model2info_running_on_matlab_engine( model, fileID, options )
% MODEL2INFO_RUNNING_ON_MATLAB_ENGINE List basic information about an SLML model.
%
% model   model structure (required)
% options (optional) options for modifying information displayed

% Ivan E. Cao-Berg, Xin Lu
%
% Copyright (C) 2019 Murphy Lab
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.NAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'name' )
    fprintf( fileID,'%%%%%% Model name\n' );
    fprintf( fileID, ['''' model.name '''' ';\n'] );
else
    fprintf( fileID,'%%%%%% Model name\n' );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.ID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'id' )
    fprintf( fileID,'%%%%%% Model ID\n' );
    fprintf( fileID, ['''' model.id '''' ';\n'] );
else
    fprintf( fileID,'%%%%%% Model ID\n' );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIMENSIONALITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'dimensionality' )
    fprintf( fileID,'%%%%%% Dimensionality\n' );
    fprintf( fileID, ['''' model.dimensionality ''''  ';\n'] );
else
    fprintf( 'Dimensionality: UNSET;\n' );
    warning( 'This is an invalid model' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DOCUMENTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'documentation' )
    fprintf( fileID,'%%%%%% Documentation\n' );
    fnames = fieldnames( model.documentation );
    fprintf( fileID,'%% Fields and values present in the documentation structure\n\n' );
    for i=1:1:length(fnames)
        fname = fnames{i};
        try
            fprintf( fileID, ['%% ' fname '\n'] );
            fprintf( fileID, ['''' eval([ 'model.documentation.' fname]) ''''  ';\n'] );
        catch err
            disp('Unable to write parameters to report');
            getReport( err );
        end
        clear fname

    end
    clear fnames
end

try
    x=model.nuclearShapeModel;
    fprintf( fileID,'\n%%%%%% Nuclear shape model information\n' );
    model2info2(model.nuclearShapeModel, fileID )
catch
    fprintf( fileID,'\n%%%%%% Nuclear shape model information\n' );
    fprintf( fileID,'%% The model file does not contain an instance of the submodel.\n' );
end

try
    x=model.cellShapeModel;
    disp(' ');
    fprintf( fileID,'\n%%%%%% Cell shape model information\n' );
    model2info2(model.cellShapeModel, fileID )
catch
    fprintf( fileID,'\n%%%%%% Cell shape model information\n' );
    fprintf( fileID,'%% The model file does not contain an instance of the submodel.\n' );
end

try
    x=model.proteinModel;
    disp(' ');
    fprintf( fileID,'\n%%%%%% Protein pattern model information\n' );
    model2info2(model.proteinModel, fileID );
catch
    fprintf( fileID,'\n%%%%%% Protein pattern model information\n' );
    fprintf( fileID,'%% The model file does not contain an instance of the submodel.\n' );
end

try
    x=model.proteinModel.spharm_obj_model.cellShapeModel;
    disp(' ');
    fprintf( fileID,'\n%%%%%% Spharm_obj model information\n' );
    model2info2(model.proteinModel.spharm_obj_model.cellShapeModel, fileID )
catch    
    fprintf( fileID,'%% The model file does not contain an instance of the submodel.\n' );
end

try
    x=model.proteinModel.ppm;
    disp(' ');
    fprintf( fileID,'\n%%%%%% ppm model information\n' );
    model2info3(model.proteinModel.ppm, fileID )
catch    
    fprintf( fileID,'%% The model file does not contain an instance of the submodel.\n' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATASET Information %%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(model, 'dataset')
	answer = datasetinfo2html(model, fileID);
	if isfield(model.dataset, 'segmentation')
		ans_seg = segmentation2html(model, fileID, 'html');
	end
end
if isfield(options, 'labels')
    labels = options.labels;
else 
    labels = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS DIFFEOMORPHIC? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_diffeomorphic( model )
    fprintf(fileID,'%%%%%% Additional Figures\n');
    f = figure('visible','off');
    nimgs = size(model.cellShapeModel.positions,1);
    options.plot_dims = [1,2];
    options.subsize = 1000;
    showShapeSpaceFigure( model, labels, options );
    saveas( f, 'show_shape_space.png', 'png' );
    I = imread( 'show_shape_space.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'show_shape_space_thumbnail.png' );

    fprintf( fileID, '%%\n' );
    fprintf( fileID, '%% <html>\n' );
    fprintf( fileID, '%% <figure>\n' );
    fprintf( fileID, '%% <a href="show_shape_space.png" ><img src="show_shape_space_thumbnail.png" /></a>\n' );
    fprintf( fileID, '%% <figcaption>Figure 1. Shape space for input model. </figcaption>\n' );
    fprintf( fileID, '%% </figure>\n' );
    fprintf( fileID, '%% </html>\n' );
    fprintf( fileID, '%%\n' );

    movefile( '*.png', './html' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS PCA FRAMEWORK? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_pca_framework( model )
    fprintf(fileID,'%%%%%% Additional Figures\n');
    f = figure('visible','off');
    showPCAShapeSpaceFigure( model, labels, options );
    saveas( f, 'show_shape_space.png', 'png' );
    I = imread( 'show_shape_space.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'show_shape_space_thumbnail.png' );

    fprintf( fileID, '%%\n' );
    fprintf( fileID, '%% <html>\n' );
    fprintf( fileID, '%% <figure>\n' );
    fprintf( fileID, '%% <a href="show_shape_space.png" ><img src="show_shape_space_thumbnail.png" /></a>\n' );
    fprintf( fileID, '%% <figcaption>Figure 1. Shape space for input model. </figcaption>\n' );
    fprintf( fileID, '%% </figure>\n' );
    fprintf( fileID, '%% </html>\n' );
    fprintf( fileID, '%%\n' );

    movefile( '*.png', './html' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS CLASSIC MODEL? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_classic_model( model )
    fprintf(fileID,'%%%%%% Additional Figures\n');
    %f = figure('visible','off');
    fprintf( fileID,'%% There are no additional figures for this model at the moment.\n' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS SPHARM MODEL? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(model,'proteinModel')
    if isfield(model.proteinModel,'spharm_obj_model')
        model=model.proteinModel.spharm_obj_model;
    end
end
if is_spharm_model( model )
    fprintf(fileID,'%%%%%% Additional Figures\n');
    showevol = true;
    if isfield(options,'shape_evolution')
        if strcmpi(options.shape_evolution,'none')
            showevol = false;
        end
    end
    if showevol
        f = figure('visible','off');
        show_SPHARM_RPDM_shape_evolution_figure(model);
        saveas( gcf, 'show_shape_evolution.png', 'png' );
        I = imread( 'show_shape_evolution.png' );
        I = imresize( I, 0.50 );
        imwrite( I, 'show_shape_evolution_thumbnail.png' );
        fprintf( fileID, '%%\n' );
        fprintf( fileID, '%% <html>\n' );
        fprintf( fileID, '%% <figure>\n' );
        fprintf( fileID, '%% <a href="show_shape_evolution.png" ><img src="show_shape_evolution_thumbnail.png" /></a>\n' );
        fprintf( fileID, '%% <figcaption>Figure 1. Shape evolution. </figcaption>\n' );
        fprintf( fileID, '%% </figure>\n' );
        close all
    end
    f = figure('visible','off');
    show_SPHARM_RPDM_Shape_Space_Figure(model.cellShapeModel,labels,options);
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

    movefile( '*.png', './html' );
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

    movefile( '*.png', './html' );
end
end

function model2info2( model, fileID )

fprintf( fileID,'%% Fields and values present in the submodel.\n\n' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.NAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'name' )
    fprintf( fileID, ['\n%% model.name\n'] );
    fprintf( fileID, ['''' model.name '''' ';\n'] );
else
    fprintf( fileID, ['%% model.name\n'] );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.ID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'id' )
    fprintf( fileID, ['%% model.id\n'] );
    fprintf( fileID, ['''' model.id '''' ';\n'] );
else
    fprintf( fileID, ['%% model.id\n'] );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'class' )
    fprintf( fileID, ['%% model.class\n'] );
    fprintf( fileID, ['''' model.class '''' ';\n'] );
else
    fprintf( fileID, ['%% model.class\n'] );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.TYPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'type' )
    fprintf( fileID, ['%% model.type\n'] );
    fprintf( fileID, ['''' model.type '''' ';\n'] );
else
    fprintf( fileID, ['%% model.type\n'] );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.RESOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'resolution' )
    fprintf( fileID, ['%% model.resolution\n'] );
    fprintf( fileID, ['''' mat2str(model.resolution) '''' ';\n'] );
else
    fprintf( fileID, ['%% model.resolution\n'] );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end
fprintf( fileID, ['\n'] );
end

function model2info3( model, fileID )

fprintf( fileID,'%% Fields and values present in the submodel.\n\n' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.NAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'cellmembrane' )
    fprintf( fileID, ['\n%% model.cellmembrane.likelihood\n'] );
    fprintf( fileID, ['''' model.cellmembrane.likelihood '''' ';\n'] );
    fprintf( fileID, ['\n%% model.cellmembrane.param\n'] );
    fprintf( fileID, ['''' model.cellmembrane.param '''' ';\n'] );
else
    fprintf( fileID, ['%% model.cellmembrane\n'] );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.ID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'nuclear' )
    fprintf( fileID, ['%% model.nuclear.likelihood\n'] );
    fprintf( fileID, ['''' model.nuclear.likelihood '''' ';\n'] );
    fprintf( fileID, ['%% model.nuclear.param\n'] );
    fprintf( fileID, ['''' model.nuclear.param '''' ';\n'] );
else
    fprintf( fileID, ['%% model.param\n'] );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield( model, 'full' )
    fprintf( fileID, ['%% model.full.likelihood\n'] );
    fprintf( fileID, ['''' model.full.likelihood '''' ';\n'] );
    fprintf( fileID, ['%% model.full.param\n'] );
    fprintf( fileID, ['''' model.full.param '''' ';\n'] );
else
    fprintf( fileID, ['%% model.full\n'] );
    fprintf( fileID, ['''UNSET''' ';\n'] );
end
fprintf( fileID, ['\n'] );
end
