function model2info_v2( model, fileID, options )
% MODEL2INFO_v2 List basic information about an SLML model.
%
% model   model structure (required)
% options (optional) options for modifying information displayed

% Ivan E. Cao-Berg, Xin Lu, Ted Zhang, Robert F. Murphy
%
% Copyright (C) 2019, 2020 Murphy Lab
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
% 2020/10/16 R. F. Murphy added printout of Hausdorff distances for
%                       SPHARM-RPDM models
% 2021/03/05 R.F. Murphy modified to use new name for hd values
% 2021/04/24 R.F. Murphy merged in Serena's show_spatial_distribution
% 2022/08/15 R.F. Murphy add hausdorff distance plots
% 2023/02/16 R.F. Murphy make figures invisible in case running deployed
% 2023/03/20 R.F. Murphy fix saving of shape space thumbnails
% 2023/04/24 R.F. Murphy add comment indicators to 3/20/2023 tag; add stats
%                           and plots for jaccard indices
% 2023/05/08 R.F. Murphy pass fileid to spharm_rpdm_sample_or_reconstruct_images_figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.NAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header2html( fileID, 'Model name');

if isfield( model, 'name' )
    text2html( fileID, model.name );
else
    text2html( fileID, 'UNSET' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.ID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header2html( fileID, 'Model ID');

if isfield( model, 'id' )
    text2html( fileID, model.id);
else
    text2html( fileID, 'UNSET');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIMENSIONALITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header2html( fileID, 'Dimensionality');

if isfield( model, 'dimensionality' )
    text2html( fileID, model.dimensionality);
else
    text2html( fileID, 'UNSET');
    warning( 'This is an invalid model' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DOCUMENTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header2html( fileID, 'Documentation');

if isfield( model, 'documentation' )
    fnames = fieldnames( model.documentation );
    text2html( fileID, 'Fields and values present in the documentation structure');
    for i=1:1:length(fnames)
        fname = fnames{i};
        try
            text2html( fileID, fname);          
            %not sure what this is doing
            text2html( fileID, [eval([ 'model.documentation.' fname])] );        
        catch err
            disp('Unable to write documentation to report');
            getReport( err );
        end
        clear fname

    end
    clear fnames
end

try
    x=model.nuclearShapeModel;
    header2html( fileID, 'Nuclear shape model information' );
    model2info2(model.nuclearShapeModel, fileID )
catch
    header2html( fileID, 'Nuclear shape model information' );
    text2html( fileID, 'The model file does not contain an instance of the submodel');
end

try
    x=model.cellShapeModel;
    header2html( fileID, 'Cell shape model information' );
    model2info2(model.cellShapeModel, fileID )
catch
    header2html( fileID, 'Cell shape model information' );
    text2html( fileID, 'The model file does not contain an instance of the submodel');
end

try
    x=model.proteinModel;
    header2html( fileID, 'Protein pattern model information' );
    model2info2(model.proteinModel, fileID );
catch
    header2html( fileID, 'Protein pattern model information' );
    text2html( fileID, 'The model file does not contain an instance of the submodel');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATASET Information %%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(model, 'dataset')
	answer = datasetinfo2html(model, fileID);
	if isfield(model.dataset, 'segmentation')
        warning('Not printing segmentations.')
		%ans_seg = segmentation2html(model, fileID, 'html');
        ans_seg = false;
	end
end

if isfield(options, 'labels')
    labels = options.labels;
else 
    labels = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS DIFFEOMORPHIC? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_diffeomorphic( model )
    
    header2html( fileID, 'Additional Figures' );
    f = figure('visible','off');
    nimgs = size(model.cellShapeModel.positions,1);
    options.plot_dims = [1,2];
    options.subsize = 1000;
    showShapeSpaceFigure( model, labels, options );
    saveas( f, 'show_shape_space.png', 'png' );
    fP = f.Position;
    f.Position = [fP(1) fP(2) round(fP(3)/2) round(fP(4)/2)];
    saveas( f, 'show_shape_space_thumbnail.png', 'png' ); %3/20/2023
    close(f)
    
    img2html(fileID, 'show_shape_space.png', 'show_shape_space_thumbnail.png', 'Diffeomorphic Model Shape space');
    movefile( '*.png', './html' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS PCA FRAMEWORK? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_pca_framework( model )
    header2html( fileID, 'Additional Figures' );
    f = figure('visible','off');
    showPCAShapeSpaceFigure( model, labels, options );
    saveas( f, 'show_shape_space.png', 'png' );
    fP = f.Position;
    f.Position = [fP(1) fP(2) round(fP(3)/2) round(fP(4)/2)];
    saveas( f, 'show_shape_space_thumbnail.png', 'png' ); %3/20/2023
    close(f)
    
    img2html(fileID, 'show_shape_space.png', 'show_shape_space_thumbnail.png', '2D PCA Model Shape space');
    movefile( '*.png', './html' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS CLASSIC MODEL? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_classic_model( model )
    % nothing to output as of now
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS SPHARM MODEL? %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
spharm_obj = false;
if(isfield(model,'proteinModel'))
    if(isfield(model.proteinModel,'spharm_obj_model'))
        model=model.proteinModel.spharm_obj_model;
        spharm_obj = true;
    end
end
if is_spharm_model( model )
    text2html( fileID, '');
    header2html( fileID, 'SPHARM-RPDM model quality');

    try
    mdl = model.cellShapeModel;
    if isfield(mdl,'hausdorff_distances') 
        if isfield(model,'cellShapeModel') 
            text2html(fileID,strcat("Cell:", hd_stats(mdl.hausdorff_distances(1,:),'Hausdorff distances')));
            if isfield(model,'nuclearShapeModel')
                text2html(fileID,strcat("Nuc:", hd_stats(mdl.hausdorff_distances(2,:),'Hausdorff distances')));
            end
        else if isfield(model,'nuclearShapeModel')
            text2html(fileID,strcat("Nuc:", hd_stats(mdl.hausdorff_distances(1,:),'Hausdorff distances')));
            end
        end
        plot_hd_or_ji(fileID,model,mdl.hausdorff_distances,'hausdorff distances','hausdorff_distances')
    end
    catch
       text2html(fileID,'Unable to display hausdorff distance statistics and plots') 
    end
    try
        if isfield(mdl,'jaccard_indices')
            if isfield(model,'cellShapeModel') 
                text2html(fileID,strcat("Cell:", hd_stats(mdl.jaccard_indices(1,:),'Jaccard indices')));
            if isfield(model,'nuclearShapeModel')
                text2html(fileID,strcat("Nuc:", hd_stats(mdl.jaccard_indices(2,:),'Jaccard indices')));
            end
        else if isfield(model,'nuclearShapeModel')
            text2html(fileID,strcat("Nuc:", hd_stats(mdl.jaccard_indices(1,:),'Jaccard indices')));
            end
        end
        plot_hd_or_ji(fileID,model,mdl.jaccard_indices,'jaccard indices','jaccard_indices')
    end
    catch
       text2html(fileID,'Unable to display jaccard index statistics and plots') 
    end
       
    
    header2html( fileID, 'Additional Figures');
    showevol = true;
    if isfield(options,'shape_evolution')
        if strcmpi(options.shape_evolution,'none')
            showevol = false;
        end
    else
        options.shape_evolution = 'first two shapes';
    end
    captionbase = 'SPHARM-RPDM Cell Shape Model: ';
    if spharm_obj captionbase = 'SPHARM-RPDM Protein Object Model: '; end
    if showevol
        spharm_rpdm_sample_or_reconstruct_images_figure(model,fileID,options);
    end
    
    f = figure('visible','off');
    show_SPHARM_RPDM_Shape_Space_Figure(model.cellShapeModel,labels,options);
    saveas( f, 'show_shape_space.png', 'png' );
    fP = f.Position;
    f.Position = [fP(1) fP(2) round(fP(3)/2) round(fP(4)/2)];
    saveas( f, 'show_shape_space_thumbnail.png', 'png' ); %3/20/2023
    img2html(fileID, 'show_shape_space.png', ...
        'show_shape_space_thumbnail.png', ...
        [captionbase 'Shape space']);
    close(f)
    
    %Spatial Distribution
    if isfield(model, 'spatial')
        show_spatial_distribution(model.spatial,fileID);
    end

    movefile( '*.png', './html' );    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IS TCell MODEL? %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_tcell_model( model )
    header2html( fileID, 'Additional Figures');
    ShowTcellEnrichment(model.model_filename, options);
    img_name=ml_ls('*model_mean_enrichment_plot_mdl*.png');
    I = imread( img_name{1} );
    imwrite( I, 'show_enrichment.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'show_enrichment_thumbnail.png' );

    img2html(fileID, 'show_enrichment.png', 'show_enrichment_thumbnail.png', 'T Cell Model Shape enrichment plot');
    movefile( '*.png', './html' );

    close all
end

end

function model2info2( model, fileID )

text2html( fileID, 'Fields and values present in the submodel');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.NAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header2html( fileID, 'Model Name');

if isfield( model, 'name' )
    text2html( fileID, model.name);
else
    text2html( fileID, 'UNSET');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.ID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header2html( fileID, 'ID');

if isfield( model, 'Model ID' )
    text2html( fileID, model.id);
else
    text2html( fileID, 'UNSET');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header2html( fileID, 'Model Class');

if isfield( model, 'class' )
    text2html( fileID, model.class);
else
    text2html( fileID, 'UNSET');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.TYPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header2html( fileID, 'Model Type');

if isfield( model, 'type' )
    text2html( fileID, model.type);
else
    text2html( fileID, 'UNSET');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.RESOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
header2html( fileID, 'Model Resolution');

if isfield( model, 'resolution' )
    text2html( fileID, mat2str(model.resolution));
else
    text2html( fileID, 'UNSET');
end

end


function model2info3( model, fileID )
text2html( fileID, 'Fields and values present in the submodel');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.Cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header2html( fileID, 'Model cell membrane');

if isfield( model, 'cellmembrane' )
    text2html( fileID, model.cellmembrane.likelihood);
    text2html( fileID, model.cellmembrane.param);
else
     text2html( fileID, 'UNSET');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.Nuclear %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header2html( fileID, 'Model Nuclear');

if isfield( model, 'nuclear' )
    text2html( fileID, model.nuclear.likelihood);
    text2html( fileID, model.nuclear.param);
else
      text2html( fileID, 'UNSET');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL.Full %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header2html( fileID, 'Model Full');

if isfield( model, 'full' )
    text2html( fileID, model.nuclear.likelihood);
    text2html( fileID, model.nuclear.param);
else
      text2html( fileID, 'UNSET');
end

end

function outstr = hd_stats(hd,tag)

[hdmin,hdmini] = min(hd);
[hdmax,hdmaxi] = max(hd);
hdavg = mean(hd);
hdsd = std(hd);
outstr = strcat(tag," avg=",num2str(hdavg), ...
    ", sd=", num2str(hdsd), ", min=", num2str(hdmin), " [", num2str(hdmini), ...
    "], max=", num2str(hdmax), " [", num2str(hdmaxi), "]");
end

function plot_hd_or_ji(fileID,model,distnces,tag,tag2)

if isfield(model,'cellShapeModel')
    f = figure('visible','off');
    histogram(log10(distnces(1,:)),25,'FaceColor','red')
    xlabel(['Log10 ' tag ' between original and SPHARM-RPDM model']);
    ylabel('Frequency')
    if isfield(model,'nuclearShapeModel')
        hold on
        histogram(log10(distnces(2,:)),25,'FaceColor','green')
        legend('cells','nuclei')
        saveas( f, [tag2 '_histogram.png'], 'png' );
        img2html(fileID, [tag2 '_histogram.png'], ...
        [tag2 '_histogram.png'], ...
        [tag ' for reconstructions']);
        hold off
        % now plot scatter of nuclear vs cell
        f = figure('visible','off');
        plot(log10(distnces(2,:)), ...
            log10(distnces(1,:)),'d')
        xlabel(['log10 ' tag ' for cell shape'])
        ylabel(['log10 ' tag ' for nuclear shape'])
        saveas( f, [tag2 '_scatter.png'], 'png' );
        img2html(fileID, [tag2 '_scatter.png'], ...
        [tag2 '_scatter.png'], ...
        tag);
    else
        legend('cells')
        saveas( f, [tag2 '_histogram.png'], 'png' );
        img2html(fileID, [tag2 '_histogram.png'], ...
        [tag2 '_histogram.png'], ...
        [tag ' for reconstructions']);
    end
elseif isfield(model,'nuclearShapeModel')
    f = figure('visible','off');
    histogram(log10(distnces(1,:)),25,'FaceColor','green')
    xlabel(['Log10 ' tag ' between original and SPHARM-RPDM model']);
    ylabel('Frequency')
    legend('nuclei')
    saveas( f, [tag2 '_histogram.png'], 'png' );
    img2html(fileID, [tag2 '_histogram.png'], ...
    [tag2 '_histogram.png'], ...
    [tag ' for reconstructions']);
end
end