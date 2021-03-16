function [ nucimg, cellimg, outres, options ] = model2framework( model, options )
% MODEL2FRAMEWORK Helper method that generates a 2D/3D framework from a
% valid SLML model.
%
% The models that are currently supported are
% (1) 2D spline model of the nucleus + 2D ratio model of the cell membrane
% (2) 3D spline model of the nucleus + 3D ratio model of the cell membrane
% (3) 2D/3D diffeomorphic model of the nucleus and cell membrane

% Author: Ivan E. Cao-Berg (icaoberg@cs.cmu.edu)
%
% Copyright (C) 2012-2019 Murphy Lab
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

% April 7, 2012 R.F. Murphy Use modified function for synthesizing cell
% shape; adding additional models for height ratios and nuclear position
%
% September 28, 2012 I. Cao-Berg Added documentation and checking of input
% arguments
%
% October 1, 2012 I. Cao-Berg Added diffeomorphic model structure so that
% the method knows what to do when a model of this type is added. I also
% added a default value so that if the model type is not recognized, the
% method returns empty nuclear and cell images
%
% October 1, 2012 I. Cao-Berg Added helper method that synthesizes an
% instance from a diffeomorphic model
%
% November 1, 2012 I. Cao-Berg Added option to parameter structure that allows
% the selection between two methods when synthesizing instances from diffeomorphic
% models. This option is meant for testing, and it is not documented in top
% level methods
%
% November 5, 2012 I. Cao-Berg Fixed bug in method that synthesizes images from
% diffeomorphic models that was returning an empty image due to wrong index
%
% November 6, 2012 I. Cao-Berg Updated synthesis code to use Euler integrator
%
% January 20, 2013 T. Buck Updated method to match new options from newly
% trained diffeomorphic model
%
% February 4, 2013 D. Sullivan made file return a resolution associated
% with the cell and nuclear shapes. If there is none specified, prints a
% warning and returns the outres field as blank.
%
% February 22, 2013 D. Sullivan made the nuclear shape synthesis dependent
%                               on the object model resolution
%
% February 24, 2013 D. Sullivan Moved resolution checks to model2img
%
% March 3, 2013 D. Sullivan Added resolution outres to the 2D synthesis
%
% March 5, 2013 D. Sullivan Added check for an existing framework
%
% May 18, 2013 I. Cao-Berg Updated method so if framework fails to
% synthesize it returns an empty framework
%
% July 1, 2013 I. Cao-Berg Updated method so that when parameter.resolution
% is not present, it will inherit the resolution of the model. This will
% only happen in the 2D case where images are not downsampled
%
% April 7, 2012 R.F. Murphy Use modified function for synthesizing cell
% shape; adding additional models for height ratios and nuclear position
%
% June 25, 2013 R.F. Murphy Added param.spherical_cell
%
% March 17, 2014 I. Cao-Berg Updated the model type for 3D vesicular
% models. It used to be 'medial axis' and it should have been 'cylindrical_surface'
%
% April 20, 2015 I. Cao-Berg Fixed bug in method that would query a field
% even if field was not present in structure.
%
% February 1, 2016 I. Cao-Berg Replaced tp_imbox calls with cropImgND
%
% March 8, 2016 X. Ruan change param to options, fix bug for debug.
%
% Feb 17, 2018 X. Ruan add code for 2D PCA model
%
% May 20, 2019, X. Ruan add resolution for SPHARM-RPDM method synthesis

if nargin > 2
    error('Wrong number of input arguments.' );
end

if ~exist( 'options', 'var' )
    options = [];
end

%icaoberg 10/1/2012
options = ml_initparam( options, struct( ...
    'debug', false, ...
    'display', false, ...
    'verbose', false ...
    ));

nucimg = [];
cellimg = [];
outres = [];
options.nucmesh = [];
options.cellmesh = [];

%icaoberg 7/1/2013
if length(model) > 1 || ( isa( model, 'cell' ) && length(model) == 1 )
    %icaoberg 9/28/2012
    if options.verbose && options.debug
        warning(['More than one model has been found. Synthesizing framework ' ...
            ' from the first model in the cell array']);
    end
    model = model{1};
end

%yajing 7/2/2013 Added field spherical_cell from bob's update
options = ml_initparam( options, struct( 'spherical_cell', false ) );

if isfield( model, 'dimensionality' )
    dimensionality = model.dimensionality;
else
    %icaoberg 9/28/2012
    disp( 'Unable to set dimensionality. Exiting method' );
    return
end

%icaoberg 7/1/2013
if ~isfield( options, 'synthesis' ) || ...
        ( ~strcmpi( options.synthesis, 'all' ) && ...
        ~strcmpi( options.synthesis, 'framework' ) && ...
        ~strcmpi( options.synthesis, 'nucleus' ) && ...
        ~strcmpi( options.synthesis, 'cell' ) )
    disp( 'Unrecognized synthesis flag. Setting to default value' );
    options.synthesis = 'all';
end

switch lower(dimensionality)
    case '2d'
        %generate 2D framework

        %generate cell framework
        %yajing 7/2/2013 Added spherical_cell case from bob's update
        if options.spherical_cell
            [nucimg,cellimg] = rm_gencellcomp2D_circle( model, options);
        elseif isfield( model.nuclearShapeModel, 'type' ) && ...
                isfield( model.cellShapeModel, 'type' ) && ...
                strcmpi( model.nuclearShapeModel.type, 'diffeomorphic' ) && ...
                strcmpi( model.cellShapeModel.type, 'diffeomorphic' )

            %use the same code for 3D, returns RGB image
            [nucimg, cellimg, options] = model2diffeomorphicInstance( model, options );
            if ~isfield(options, 'resolution')
                options.resolution = model.proteinModel.resolution;
                outres = model.proteinModel.resolution;
            end
            nucEdge = bwboundaries( double(sum(nucimg,3)) > 0);
            cellEdge = bwboundaries(double(sum(nucimg,3)+sum(cellimg,3)) > 0);
            cellEdge = cellEdge{1};
            nucEdge = nucEdge{1};
            nucimg = zeros(size(nucimg, 1), size(nucimg, 2));
            cellimg = zeros(size(cellimg, 1), size(cellimg, 2));
            for i = 1 : size( nucEdge, 1)
                nucimg(nucEdge(i, 1), nucEdge(i, 2)) = 1;
            end
            for i = 1 : size( cellEdge, 1)
                cellimg(cellEdge(i, 1), cellEdge(i, 2)) = 1;
            end
            % xruan 02/17/2018 add synthesize method for PCA model
        elseif  isfield( model.nuclearShapeModel, 'type' ) && ...
                isfield( model.cellShapeModel, 'type' ) && ...
                strcmpi( model.nuclearShapeModel.type, 'pca' ) && ...
                strcmpi( model.cellShapeModel.type, 'pca' )
            if isfield( options.model.pca, 'pca_synthesis_method' )
                if strcmp(options.model.pca.pca_synthesis_method, 'reconstruction')
                    [nucimg, cellimg] = pca_reconstruct_training_images( model, options );
                elseif strcmp(options.model.pca.pca_synthesis_method, 'random_sampling')
                    [nucimg, cellimg] = pca_random_sample_images( model, options );
                end
            else
                disp('Synthesis method not set. Setting to random_sampling')
                [nucimg, cellimg] = pca_random_sample_images( model, options );
            end
        else
            [nucimg,cellimg] = ml_gencellcomp2D( model, options );
        end

        %icaobeg 5/18/2013
        %make sure the nuclear edge is not empty
        %the nuclear edge can never be empty
        if ( strcmpi( options.synthesis, 'nucleus' ) || ...
                strcmpi( options.synthesis, 'framework' ) || ...
                strcmpi( options.synthesis, 'all' ) ...
                ) && isempty( nucimg )
            disp( 'Nuclear image is empty. Returning empty framework.' );
            outres = [];
            return
        end

        if ( strcmpi( options.synthesis, 'cell' ) || ...
                strcmpi( options.synthesis, 'framework' ) || ...
                strcmpi( options.synthesis, 'all' ) ...
                ) && isempty( cellimg )
            disp( 'Cell image is empty. Returning empty framework.' );
            outres = [];
            return
        end

        if options.framework_cropping
            if strcmpi( options.synthesis, 'framework' ) || ...
                    strcmpi( options.synthesis, 'all' )
                [ croppedcellimg, cropbounds ] = cropImg( cellimg );
                nucimg = nucimg(cropbounds(1):cropbounds(2),cropbounds(3):cropbounds(4),:);
                cellimg = cellimg(cropbounds(1):cropbounds(2),cropbounds(3):cropbounds(4),:);
            else
                %icaoberg 02/01/2016
                [ croppednucimg, cropbounds ] = cropImgND( nucimg );
                nucimg = nucimg(cropbounds(1):cropbounds(2),cropbounds(3):cropbounds(4),:);
            end

            if isstruct(options.nucmesh)
                options.nucmesh.vertices(:, 1) = options.nucmesh.vertices(:, 1) - (cropbounds(3)-1);
                options.nucmesh.vertices(:, 2) = options.nucmesh.vertices(:, 2) - (cropbounds(1)-1);
                % options.nucmesh.vertices(:, 3) = options.nucmesh.vertices(:, 3) - cropbounds(5);
            end
            if isstruct(options.cellmesh)
                options.cellmesh.vertices(:, 1) = options.cellmesh.vertices(:, 1) - (cropbounds(3)-1);
                options.cellmesh.vertices(:, 2) = options.cellmesh.vertices(:, 2) - (cropbounds(1)-1);
                % options.cellmesh.vertices(:, 3) = options.cellmesh.vertices(:, 3) - cropbounds(5);
            end
        else
            if strcmpi( options.synthesis, 'framework' ) || ...
                    strcmpi( options.synthesis, 'all' )
                cropbounds = [1, size(cellimg, 1), 1, size(cellimg, 2), 1, size(cellimg, 3)];
            else
                cropbounds = [1, size(nucimg, 1), 1, size(nucimg, 2), 1, size(nucimg, 3)];
            end
        end
        options.framework_xrange = cropbounds(3:4);
        options.framework_yrange = cropbounds(1:2);
        % options.framework_zrange = cropbounds(5:6);

        %D. Sullivan 3/3/12
        %Set the output resolution
        if isfield( options,'resolution' )
            if isfield( options.resolution,'objects' )
                outres = options.resolution.objects;
            end
            outres = options.resolution;
        else
            %icaoberg 7/1/2013
            disp( 'Generated framework from 2D models are not downsampled. Inheriting resolution from model.' );
            outres = model.nuclearShapeModel.resolution;
        end
    case '3d'
        %D. Sullivan 3/5/13
        disp('Check if a framework is present (as in from raw data)');
        if ( isfield( options, 'cell' ) && isfield( options.cell, 'instance' ) ) && ...
                (isfield( options, 'nucleus' ) && isfield( options.nucleus, 'instance' ))
            if all(size(options.cell.instance) == size(options.nucleus.instance))
                disp('Setting nuclear membrane instance from raw data');
                nucimg = options.nucleus.instance;

                disp('Setting cell membrane instance from raw data');
                if isfield( options.cell, 'instance' )
                    cellimg = options.cell.instance;
                end

                disp('Setting framework resolution');
                if isfield( options, 'resolution' )
                    if isfield(options.resolution, 'cell' )
                        outres = options.resolution.cell;
                    end
                else
                    outres = [];
                end
                return
            end

            if options.debug
                warning('Cell and nucleus found, but the images were not the same size, synthesizing new framework');
            end
        end

        %icaoberg march 17, 2014
        if strcmpi( options.synthesis, 'nucleus' ) && ...
                strcmpi( model.nuclearShapeModel.type, 'cylindrical_surface' )
            %generate 3D framework

            %icaoberg 8/7/2013
            %resolution is always needed
            options.resolution.nucleus = model.nuclearShapeModel.resolution;

            %D. Sullivan 2/22/13 added param structure to pass the cell and
            %object resolutions
            instance = tp_genspsurf(model.nuclearShapeModel,options);
            %instance = tp_genspsurf(model.nuclearShapeModel);

            %disp('Generating nuclear shape instance');
            %Yajing Tang 7/2/2013 Added spherical_cell case from Bob's update
            tp_gennucshape_result = tp_gennucshape(instance,options);
            nucimg = tp_gennucshape_result.nucimg;
            nucsurf = tp_gennucshape_result.nucsurf;
            options.nucmesh = tp_gennucshape_result.nucmesh;
            outres = tp_gennucshape_result.outres;

            %icaoberg 02/01/2016
            if options.framework_cropping
                [ croppednucimg, cropbounds ] = cropImgND( nucimg );
                nucimg = nucimg(cropbounds(1):cropbounds(2),cropbounds(3):cropbounds(4),:);
                cellimg = [];

                if isstruct(options.nucmesh)
                    options.nucmesh.vertices(:, 1) = options.nucmesh.vertices(:, 1) - (cropbounds(3)-1);
                    options.nucmesh.vertices(:, 2) = options.nucmesh.vertices(:, 2) - (cropbounds(1)-1);
                    % options.nucmesh.vertices(:, 3) = options.nucmesh.vertices(:, 3) - cropbounds(5);
                    nucmesh = options.nucmesh;
                end
            else
                cropbounds = [1, size(nucimg, 1), 1, size(nucimg, 2), 1, size(nucimg, 3)];
            end
            options.framework_xrange = cropbounds(3:4);
            options.framework_yrange = cropbounds(1:2);
            options.framework_zrange = cropbounds(5:6);


        elseif ( strcmpi( options.synthesis, 'cell' ) || ...
                strcmpi( options.synthesis, 'framework' ) || ...
                strcmpi( options.synthesis, 'all' ) ) && ...
                strcmpi( model.nuclearShapeModel.type, 'cylindrical_surface' ) && ...
                strcmpi( model.cellShapeModel.type, 'ratio' )

            if (isfield( options, 'nucleus' )  && isfield( options.nucleus, 'instance' ))
                if options.verbose
                    disp('Generating nuclear shape from nuclear image');
                end
                [nucimg, nucsurf, outres] = generate_nuclear_shape_from_image( options.nucleus.instance, options );
                nucleus.nucimgsize = size(nucimg);
                nucleus.nucsurf = nucsurf;
                nucleus.nucimg = nucimg;
                nucleus.instance = instance;

                if options.spherical_cell
                    [cellimg,cellsurf] = rm_gensphere(1.0);
                else
                    ml_gencellshape3d_result = ml_gencellshape3d( ...
                        model.cellShapeModel, nucleus,options );
                    cellimg = ml_gencellshape3d_result.cellimg;
                    nucimg = ml_gencellshape3d_result.nucimg;
                    options.cellmesh = ml_gencellshape3d_result.cellmesh;
                end
            else
                if options.verbose
                    disp( 'Generating nuclear shape' );
                end
                instance = tp_genspsurf(model.nuclearShapeModel,options);

                %disp('Generating nuclear shape instance');
                %Yajing Tang 7/2/2013 Added spherical_cell case from Bob's update
                if options.spherical_cell
                    % this doesn't work yet
                    [nucimg,nucsurf] = rm_gensphere([1024*1.25, 1024*1.25, 5*instance.height],[0.33, 0.33, 0.33]);
                else
                    tp_gennucshape_result = tp_gennucshape(instance,options);
                    nucimg = tp_gennucshape_result.nucimg;
                    nucsurf = tp_gennucshape_result.nucsurf;
                    options.nucmesh = tp_gennucshape_result.nucmesh;
                    outres = tp_gennucshape_result.outres;
                end

                nucleus.nucimgsize = size(nucimg);
                nucleus.nucsurf = nucsurf;
                nucleus.nucimg = nucimg;
                nucleus.instance = instance;

                if options.verbose
                    disp('Generating cell shape');
                end
                if options.spherical_cell
                    [cellimg,cellsurf] = rm_gensphere(1.0);
                else
                    ml_gencellshape3d_result = ml_gencellshape3d( ...
                        model.cellShapeModel, nucleus,options );
                    cellimg = ml_gencellshape3d_result.cellimg;
                    nucimg = ml_gencellshape3d_result.nucimg;
                    options.cellmesh = ml_gencellshape3d_result.cellmesh;
                end
            end
        elseif strcmpi( model.nuclearShapeModel.type, 'diffeomorphic' ) && ...
                strcmpi( model.cellShapeModel.type, 'diffeomorphic' )
            %icaoberg 10/1/2012

            outres = options.resolution.cell;
            [nucimg, cellimg, options] = model2diffeomorphicInstance( model, options );

        elseif strcmpi( model.nuclearShapeModel.class, 'csgo' ) && ...
                strcmpi( model.cellShapeModel.class, 'csgo' ) && ...
                strcmpi( model.nuclearShapeModel.type, 'half_ellipsoid' ) && ...
                strcmpi( model.cellShapeModel.type, 'half_ellipsoid' )

            outres = [];
            disp( 'Generating half-ellipsoid geometry on model' );
            [nucimg, cellimg] = ellipsoid_geometry( model );

            % xruan 09/17/2018 add synthesize method for SPHARM-RPDM model
        elseif  isfield( model.nuclearShapeModel, 'type' ) && ...
                isfield( model.cellShapeModel, 'type' ) && ...
                strcmpi( model.nuclearShapeModel.type, 'spharm_rpdm' ) && ...
                strcmpi( model.cellShapeModel.type, 'spharm_rpdm' )
            if strcmp(options.spharm_rpdm.synthesis_method, 'reconstruction') || strcmp(options.spharm_rpdm.synthesis_method, 'random_sampling')
                % xruan 05/20/2019 check if there is resolution for
                % protein, if so, use 'cell' as the resolution for framework
                if strcmpi(options.synthesis, 'all') && isfield(options.resolution, 'cell')
                    options.spharm_rpdm.synthesis_resolution = options.resolution.cell;
                end

                spharm_rpdm_sample_or_reconstruct_images_result = spharm_rpdm_sample_or_reconstruct_images( model, options );
                nucimg = spharm_rpdm_sample_or_reconstruct_images_result.nucimg;
                cellimg = spharm_rpdm_sample_or_reconstruct_images_result.cellimg;
                options.nucmesh = spharm_rpdm_sample_or_reconstruct_images_result.nucmesh;
                options.cellmesh = spharm_rpdm_sample_or_reconstruct_images_result.cellmesh;
                options.spharm_rpdm.shape_space_coords = spharm_rpdm_sample_or_reconstruct_images_result.shape_space_coords;
                options.spharm_rpdm.reconst_spharm_descriptors = spharm_rpdm_sample_or_reconstruct_images_result.reconst_spharm_descriptors;
            end
            outres = options.resolution.cell;
        else
            %icaoberg 10/1/2012
            warning( 'CellOrganizer: Unrecognized model type or combination of model types.' );
            nucimg = [];
        end

        if strcmpi( options.synthesis, 'nucleus' )
            cellimg = [];
        elseif strcmpi( options.synthesis, 'cell' )
            nucimg = [];
        end

        cropbounds = [];
        if options.framework_cropping
            if ~isempty(nucimg) && ~isempty(cellimg)
                [~, cropbounds] = cropImg(nucimg+cellimg);
            elseif ~isempty(nucimg)
                [~, cropbounds] = cropImg(nucimg);
            elseif ~isempty(cellimg)
                [~, cropbounds] = cropImg(cellimg);
            end
        end

        if options.framework_cropping && ~isempty(cropbounds)
            if ~isempty(nucimg)
                nucimg = nucimg(cropbounds(1):cropbounds(2), cropbounds(3):cropbounds(4),:);
            end
            if ~isempty(cellimg)
                cellimg = cellimg(cropbounds(1):cropbounds(2), cropbounds(3):cropbounds(4),:);
            end
            if isstruct(options.nucmesh)
                options.nucmesh.vertices(:, 1) = options.nucmesh.vertices(:, 1) - (cropbounds(3)-1);
                options.nucmesh.vertices(:, 2) = options.nucmesh.vertices(:, 2) - (cropbounds(1)-1);
                % options.nucmesh.vertices(:, 3) = options.nucmesh.vertices(:, 3) - cropbounds(5);
            end
            if isstruct(options.cellmesh)
                options.cellmesh.vertices(:, 1) = options.cellmesh.vertices(:, 1) - (cropbounds(3)-1);
                options.cellmesh.vertices(:, 2) = options.cellmesh.vertices(:, 2) - (cropbounds(1)-1);
                % options.cellmesh.vertices(:, 3) = options.cellmesh.vertices(:, 3) - cropbounds(5);
            end
        else
            cropbounds = [1, size(nucimg, 1), 1, size(nucimg, 2), 1, size(nucimg, 3)];
        end
        options.framework_xrange = cropbounds(3:4);
        options.framework_yrange = cropbounds(1:2);
        options.framework_zrange = cropbounds(5:6);

    otherwise
        return
end
end
