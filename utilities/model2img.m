function [imgs, param] = model2img( models, param )
% MODEL2IMG Synthesizes a multicolor image from a set of SLML Level 1.0
% Version 1.* instances. Saves the results to an intermedia folder and returns
% the synthesized image.
%
% List Of Input Arguments Descriptions
% ----------------------- ------------
% models                  A cell array of valid and parsed SLML instances.
% param                   A structure holding the function options
%
% The shape of param is described

% List Of Parameters        Descriptions
% ------------------        ------------
% targetDirectory           (optional) Directory where the images are going to be saved. Default is current directory.
% prefix                    (optional) Filename prefix for the synthesized images. Default is 'demo'
% nimgs                     (optional) Number of synthesized images. Default is 1.
% compression               (optional) Compression of tiff, i.e. 'none', 'lzw' and 'packbits'
% verbose                   (optional) Print the intermediate steps to screen. Default is true.
% microscope                (optional) Microscope model from which we select a point spread function. Default is 'none'
% sampling.method           (optional) Can be 'disc', 'trimmed' or 'sampled'. Default is trimmed
% sampling.density          (optional) An integer. Default is empty.
% protein.cytonuclearflag   (optional) Can 'cyto', 'nucleus' or 'all'. Default is all.
%
% Example
% -------
% instances = { 'model' };
% param.targetDirectory = pwd;
% param.prefix = 'demo'
% param.nimgs = 10;
% param.compression = 'lzw';
% param.microscope = 'svi';
% param.verbose = true;
%
% >> slml2img( instances, param );

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2008-2018 Murphy Lab
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

% March 7, 2012 R. F. Murphy Change scaling to 0-255
%
% March 8, 2012 @icaoberg Documented method
%
% March 9, 2012 @icaoberg Commented out sampling density,
% changed disp to fprintf and improved logs
%
% March 14, 2012 @icaoberg Added protein.location parameter
%
% March 14, 2012 @icaoberg Added contrast stretching for every pattern
%
% March 21, 2012 @icaoberg Changed protein.location to protein.cytonuclear flag
%
% April 11, 2012 @icaoberg Parameter structure is set to empty by default
% when the number of input arguments is one. Added try/catch clause to
% determine verbose flag. Default is false.
%
% April 12, 2012 R.F. Murphy Fix calculation of nucleardist and celldist;
% fix name of microtubule model
%
% August 4, 2012 D. Sullivan Used the param.nucleus and param.cell instead
% of creating separate "filled" images to save memory.
%
% August 6, 2012 @icaoberg Updated access to temporary folder so that it access the
% correct folder if more than one temporary folder is present in the current workspace path
%
% October 1, 2012 @icaoberg Updated code so that if model2framework fails
% to make a framework and returns empty arrays, then this method returns an
% empty cell array as well as throwing a warning if the debug flag is true
%
% October 9, 2012 @icaoberg Added synthesis option that allows user to
% synthesize and return a nuclear instance ('nucleus'), a cell instance
% ('cell'), a framework instance ('framework') or make an instance with all
% the patterns in the existing model files ('all'). Default is all
%
% November 14, 2012 D. Sullivan added param.resolution.cell and
% param.resolution.objects and code to adjust frameworks to proper size for
% object synthesis.
%
% November 16, 2012 @icaoberg Updated code so that parameter structure
% can be passed to model2framework to allow model2diffeomorphic to inherit
% parameter values from slml2img
%
% January 21, 2013 D. Sullivan updated resolution framework s.t. user may
% now specify multiple object model resolutions and the output will be in
% the form of the lowest resolution (highest numbers since resolutions =
% microns/pixel)
%
% January 22, 2013 @icaoberg modified if statement so that it checks
% whether field exists before querying it
%
% January 22, 2013 @icaoberg extracts resolution information from first
% model in list of models
%
% February 20, 2013 D. Sullivan fixed resolution extraction
%
% April 29, 2013 D. Sullivan Removed duplicate code in 2D case
%
% May 8, 2013 Improved use of the param.synthesis flag for the 2D case
%
% May 13, 2013 Fixed bug where it was not setting the param.synthesis to
% the default value
%
% May 29, 2013 G. Johnson Fixed bug in 2D that broke when no protein model was
% specified
%
% July 23, 2013 D. Sullivan Added checks for param.randomwalk which
% indicates if synthesis of a walk through shape space should be done
%
% Jul 26, 2013 G. Johnson If no synthesis param is specified, synthesize only what
% is available
%
% April 28, 2014 @icaoberg Changed saving of temp results to
% use 'v7.3' flag. Since we only support versions above 7.3,
% this should not be a problem
%
% April 28, 2014 @icaoberg Bounds image to a box after
% processing to prevent out of memory error
%
% March 3, 2016 @icaoberg Included check of the param field 'synthesis'
%
% May 30, 2016 I . Cao-Berg Fixed reference to model class and type
%
% January 13, 2018 @icaoberg Cleaned up method
%
% June 11, 2018 @icaoberg Fixed issue with param.resolution.objects

%step0: check input parameters
if nargin > 2
    error( 'CellOrganizer: Wrong number of input arguments.' );
end

if nargin == 1
    param = [];
end

param = get_cellorganizer_default_parameters( 'synthesis', param );

imgs = {};
param.cellmesh = [];
param.nucmesh = [];

%icaoberg 10/29/2013
if isempty( models )
    warning( 'Input parameters models cannot be empty. Exiting method' );
    return
end

%these are the default baby parameters for cellorganizer
param = ml_initparam( param, struct( ...
    'debug', false, ...
    'display', false , ...
    'verbose', true, ...
    'fileID', [], ...
    'logger', false));

fileID = param.fileID;
logger = param.logger;

try
    sampling = param.sampling.method;
catch err
    param.sampling.method = 'disc';
end

N = [];

try
    logger = param.logger;
catch err
    logger = false;
end

%icaoberg 21/4/2014
%updated code to determine the synthesis parameter correcly
%TODO: redesign param so it doesn't overwrite
try
    synthesis = param.synthesis;
    if ~isa( synthesis, 'char' )
        if  isfield( param.synthesis, 'synthesis' )
            synthesis = param.synthesis.synthesis;
        else
            synthesis = 'all';
        end
    else
        if ~strcmpi( synthesis, 'nucleus' ) && ~strcmpi( synthesis, 'cell' ) && ...
                ~strcmpi( synthesis, 'framework' ) && ~strcmpi( synthesis, 'all' )
            if  isfield( param.synthesis, 'synthesis' )
                synthesis = param.synthesis.synthesis;
            else
                synthesis = 'all';
            end
        else
            synthesis = lower( synthesis );
        end
    end
catch
    %icaoberg 5/13/2013
    %grj 7/26/13 - if no synthesis param is specified, synthesize only what
    %is available
    if isfield(models{1}, 'proteinModel')
        param.synthesis = 'all';
    elseif isfield(models{1}, 'cellShapeModel')
        param.synthesis = 'framework';
    else
        param.synthesis = 'nucleus';
    end
    synthesis = param.synthesis;
end

%checking the validity of the input models

%initialize parameters
if isempty( models )
    disp( 'Input argument models cannot be empty' );
    imgs = [];
    return;
else
    if ~isempty( fileID )
        fprintf( fileID, '%s\n', ['Setting model dimensionality to ' models{1}.dimensionality] );
    end

    if param.verbose
        fprintf( 1, '%s\n', ['Setting model dimensionality to ' models{1}.dimensionality] );
    end

    try
        dimensionality = models{1}.dimensionality;
    catch err
        warning('Unable to set model dimensionality');
        return
    end
end

disp( 'Checking all models have the same dimensionality' );
if ~exist( param.temporary_results, 'dir' )
    mkdir( param.temporary_results );
end

for i=1:1:length(models)
    %icaoberg 10/1/2012
    if ~strcmpi( models{i}.dimensionality, dimensionality )
        %icaoberg 7/1/2013
        disp('Models have different dimensionality. Unable to synthesize image. Exiting method.');

        if ~isempty( fileID )
            fprintf( fileID, '%s\n', 'Models have different dimensionality. Unable to synthesize image.' );
        end

        return;
    end
end

switch lower(dimensionality)
    %synthesizes 2D images
    case '2d'
        if ~exist( 'param', 'var' )
            param  = [];
        end

        imgs = [];
        %the size of the images is fixed to be 1024x1024
        %D. Sullivan 4/29/13, this is redundent with the
        %ml_initparam call on ~line110 of model2framework removed the
        %model2framework one

        if ~isfield( param, 'image_size' )
            param = ml_initparam(param,struct('imageSize',[1024 1024]));
        else
            param = ml_initparam(param,struct('imageSize', param.image_size));
            param = rmfield( param, 'image_size' );
        end

        param = ml_initparam(param, ...
            struct('gentex',0,'loc','all', 'nimgs', 1));

        imgs = {};
        for i = 1:param.nimgs

            [nucEdge,cellEdge,outres,param] = model2framework( models{1}, param );

            param.resolution.cell = outres;

            %icaoberg 13/5/2012
            if strcmpi( param.synthesis, 'framework' ) || ...
                    strcmpi( param.synthesis, 'all' )

                %if you asked cellorganizer to synthesize a framework and
                %either is empty, then return with a warning
                if isempty( nucEdge ) || isempty( cellEdge )
                    if param.debug
                        warning('Unable to synthesize 2D framework.');
                        return
                    end
                end
            else
                %if you asked cellorganizer to synthesize a nucleus and it is
                %empty, then return with a warning
                if isempty( nucEdge )
                    warning('Unable to synthesize 2D nuclear shape.');
                    return
                end
            end

            %icaoberg 5/13/2013
            dilate = false;
            if dilate
                %post-process each channel
                se=strel('disk',4,4);
                if ~exist('nuctex','var')
                    nucimage = imdilate(nucEdge,se);
                else
                    nucimage = nucEdge;
                end

                %icaoberg 5/13/2013
                if strcmpi( param.synthesis, 'framework' ) || ...
                        strcmpi( param.synthesis, 'all' )
                    cellimage = imdilate(cellEdge,se);
                end
            else
                %post-process each channel
                nucimage = nucEdge;

                if strcmpi( param.synthesis, 'framework' ) || ...
                        strcmpi( param.synthesis, 'all' )
                    cellimage = cellEdge;
                end
            end

            imgs = [];
            img(:,:,1) = ml_bcimg(double(nucimage),[],[0 255]);

            if exist('cellimage', 'var') && ~isempty(cellimage) && ...
                    ( strcmpi(synthesis, 'all') || ...
                    strcmpi(synthesis, 'framework' ) )
                img(:,:,2) = ml_bcimg(double(cellimage),[],[0 255]);
            end

            if strcmpi(synthesis,'all')
                %G. Johnson 5/29/13 added check to see if protein model exists
                c = 1;
                protimage = {};
                for j=1:1:length(models)
                    %If the protein model is empty the only field will be
                    %'resolution'
                    % xruan 01/07/2016
                    % the above is not true currently, it contains
                    % resolution, dimension even if the model is empty
                    % so I changed to 3
                    % if isfield(models{j},'proteinModel') && length(fieldnames(models{j}.proteinModel)) > 1
                    if isfield(models{j},'proteinModel') && length(fieldnames(models{j}.proteinModel)) > 2
                        protimage{c} = ...
                            ml_genprotimg2D( models{j}.proteinModel, ...
                            nucEdge, cellEdge, param );
                        c = c+1;
                    end
                end

                for j=1:1:length(protimage)
                    img(:,:,2+j) = ml_bcimg(double(protimage{j}),[],[0 255]);
                end
            end
            imgs{i} = img;
        end

        save( [param.temporary_results filesep 'image.mat'], 'imgs' );
    case '3d'
        %synthesize 3D images
        %only do this check if synthesis is set to all'
        disp(['Checking synthesis flag: ' synthesis ] );

        %icaoberg 7/8/2013
        %this will verify the first model given the user only as
        if strcmpi( synthesis, 'nucleus' )
            disp( 'Checking parameters for the nuclear model from the first model in list');
        end


        %icaoberg 7/8/2013
        %we only care about checking the parameters for the framework iff
        %the synthesis flag is either set to framework or to all
        if strcmpi( synthesis, 'cell' ) || strcmpi( synthesis, 'framework' ) || strcmpi( synthesis, 'all' )
            disp( 'Checking parameters for the cell model from the first model in list');
            %the cell framework, that is the nuclear and cell shape, are
            %generated using the first model found in the list
            fprintf( 1, '%s\n', 'Generating cell framework' );
            %icaoberg 11/16/2012

            %D. Sullivan 2/24/13 moved resolution setup to above
            %Note, setting initial framework resolution off first object
            %model
            %model2framework
            %Framework(cell)
            if ~isfield( models{1}, 'cellShapeModel' )
                warning( 'First model in list did not contain a cell model. Exiting method' );
                return
            end

            if isfield(models{1}.cellShapeModel,'resolution')
                %D. Sullivan 3/6/13 check if a framework instance was
                %passed and if so use that resolution. If none was
                %specified, return an error and exit
                    param.resolution.cell = models{1}.cellShapeModel.resolution;
            else
                error('No cell resolution specified. Unable to synthesize image without this information');
            end
        end

        %icaoberg 7/8/2013
        %we only care about the resolution of objects iff the user wants to
        %synthesize objects
        if strcmpi( synthesis, 'all' )
            disp( 'Checking parameters for the protein model from the first model in list');

            %objects
            for i = 1:length(models)
                if isfield(models{i}, 'proteinModel')
                    if isfield(models{i}.proteinModel,'resolution')
                        param.resolution.objects = models{i}.proteinModel.resolution;
                    elseif isfield(param.resolution.objects)

                    else
                        error('No object resolution specified. Unable to synthesize image without this information.');
                    end
                end
            end
        end

        %icaoberg 8/7/2013
        %this will verify the models have all the right information given
        %the fact the user asked to synthesize protein patterns
        if strcmpi( synthesis, 'all' )
            disp( 'Checking all the protein models from the list');

            disp( 'Calculating number of vesicle models' );
            vesicle = {};
            for j=1:1:length(models)
                if strcmpi( models{i}.proteinModel.type, 'vesicle' )
                    vesicle{length(vesicle)+1} = j;
                end
            end

            if ~isempty( vesicle )
                fprintf( 1, '%s\n', ['At least one vesicle model was found, looking for a cell shape ' ...
                    ' and a nuclear shape model.']);
                if isfield( models{1}, 'nuclearShapeModel' ) && isfield( models{1}, 'cellShapeModel' )
                    fprintf( 1, '%s\n', 'A nuclear shape and cell shape was found in the first model.');
                end
            end

            microtubules = {};
            for j=1:1:length(models)
                if isfield(models{j}, 'proteinModel')
                    if strcmpi( models{j}.proteinModel.type, 'microtubule_growth' ) && ...
                            strcmpi( models{j}.proteinModel.class, 'network' )
                        microtubules{length(microtubules)+1} = j;
                    end
                end
            end

            if ~isempty( microtubules )
                fprintf( 1, '%s\n', ['At least one microtubule model was found, ' ...
                    'looking for a centrosome model.'] );
                centrosome = {};
                for j=1:1:length(models)
                    if strcmpi( models{j}.proteinModel.type, 'gmm' ) && ...
                            strcmpi( models{j}.proteinModel.class, 'centrosome' )
                        centrosome{length(centrosome)+1} = j;
                    end
                end
            end

            if ~isempty( microtubules ) && isempty( centrosome )
                if ~isempty( fileID )
                    fprintf( fileID, '%s\n', ['At least one microtubule model was found ' ...
                        'in your list but no centrosome models were found.'] );
                    fprintf( fileID, '%s\n', 'Closing log file' );
                    fclose( fileID );
                end

                error( [num2str(length(microtubules)) ' microtubule models where found ' ...
                    'in your list but no centrosome models found.'] ) %uncool
            end

            if exist( 'centrosome', 'var' )
                fprintf( 1, '%s\n', ['Centrosome models found:' num2str(length(centrosome))] );
            end
        end

        %at this point, we have verified the contents of the list of models
        [param.nucleus, param.cell, outres, param.synthparam] = model2framework( models{1}, param );
        param.resolution.cell = outres;
        param.instance.cell = param.cell;
        param.instance.nucleus = param.nucleus;
        
        param.cellmesh = param.synthparam.cellmesh;
        param.nucmesh = param.synthparam.nucmesh;
        
        if isfield(param.synthparam, 'spharm_rpdm')
            param.spharm_rpdm.shape_space_coords = param.synthparam.spharm_rpdm.shape_space_coords;
            param.spharm_rpdm.reconst_spharm_descriptors = param.synthparam.spharm_rpdm.reconst_spharm_descriptors;
        end

        %icaoberg 7/8/2013
        if strcmpi( param.synthesis, 'nucleus' )
            if isempty( param.nucleus )
                if param.debug
                    warning( 'Unable to synthesize 3D nucleus' );
                end

                imgs = {};
                return
            end
        elseif strcmpi( param.synthesis, 'cell' )
            if isempty( param.cell )
                if param.debug
                    warning( 'Unable to synthesize 3D cell' );
                end

                imgs = {};
                return
            end
        elseif strcmpi( param.synthesis, 'framework' )
            if isempty( param.nucleus ) || isempty( param.cell )
                if param.debug
                    warning( 'Unable to synthesize 3D framework' );
                end

                imgs = {};
                return
            end
        end

        fprintf( 1, '%s\n', 'Filling cell and nuclear images' );

        %D. Sullivan 11/14/12
        %moved resolution code from the model2instance in case the user
        %trained all object models at a single resolution
        %**this also allows the distance images to be calculated after
        %the resolution adjustment**
        %Check if the protein models were learned at a different
        %resolution than the cell shape models
        %D. Sullivan 1/21/13 moved above cell/nuc image save so the
        %images are saved at the correct resolution

        if strcmpi( param.synthesis, 'all' )
            try
                %D. Sullivan 2/20/13
                %loop through and get resolution for each object model
                param.resolution.objects = zeros(length(models),3);
                for i = 1:length(models)
                    if isfield(models{i}, 'proteinModel') && isfield(models{i}.proteinModel,'resolution')
                        param.resolution.objects(i,:) = models{i}.proteinModel.resolution;
                    elseif isfield(param.resolution,'objects')

                    else
                        error('No object resolutions specified!');
                    end
                end

                if( size( param.resolution.objects, 1 ) == 1 )
                    %D. Sullivan 7/18/13 added a function to resize so this can
                    %be used in other places
                    %Resize nucleus
                    param.nucleus = AdjustResolutions(param.nucleus,...
                        param.resolution.cell,param.resolution.objects);

                    %Resize cell
                    param.cell = AdjustResolutions(param.cell,...
                        param.resolution.cell,param.resolution.objects);

                    temp = param.resolution.cell ./ param.resolution.objects;
                    if isstruct(param.nucmesh)
                        param.nucmesh.vertices = param.nucmesh.vertices .* repmat(temp, size(param.nucmesh.vertices, 1), 1);
                    end
                    if isstruct(param.cellmesh)
                        param.cellmesh.vertices = param.cellmesh.vertices .* repmat(temp, size(param.cellmesh.vertices, 1), 1);
                    end

                    %fix the param.resolution.cell to be the new resolution to
                    %keep track of it correctly
                    param.resolution.cell = param.resolution.objects;
                    %                 param.cell =
                    %set output resolution
                    if ~isfield(param, 'outputres')
                        param.outputres = min(param.resolution.objects,1);
                    end
                    %                 %since the cell is built conditionally on the nucleus, they should
                    %                 %ALWAYS have the same resolution
                    %                 initres = param.resolution.cell;%e.g. [0.05,0.05,0.35]
                    %                 finalres = param.resolution.objects;
                    %                 param.outputres=finalres;
                    %
                    %                 %nucleus
                    %                 finalsize_x = floor(initres(1)./finalres(1).*size(param.nucleus,1));
                    %                 finalsize_y = floor(initres(2)./finalres(2).*size(param.nucleus,2));
                    %                 finalsize_z = floor(initres(3)./finalres(3).*size(param.nucleus,3));
                    %
                    %                 param.nucleus = imresize(param.nucleus,[finalsize_x finalsize_y],'bilinear');
                    %                 %need to resize the z
                    %                 param.nucleus = tp_stretch3d(param.nucleus,finalsize_z);
                    %                 param.nucleus = param.nucleus>0;
                    %
                    %                 %cell
                    %                 %Note: the recalculation of final size should be unnecessary since the
                    %                 %cell and nucleus images should always be the same size, but since the
                    %                 %arithmatic is trivially fast I recalculate to deal with potential
                    %                 %weird situations in the future(e.g. if we need the nuclear image to be
                    %                 %a smaller object that we add in to the cell image for space)DPS
                    %                 finalsize_x = floor(initres(1)./finalres(1).*size(param.cell,1));
                    %                 finalsize_y = floor(initres(2)./finalres(2).*size(param.cell,2));
                    %                 finalsize_z = floor(initres(3)./finalres(3).*size(param.cell,3));
                    %                 param.cell = imresize(param.cell,[finalsize_x finalsize_y],'bilinear');
                    %                 %need to resize the z
                    %                 param.cell = tp_stretch3d(param.cell,finalsize_z);
                    %                 param.cell = param.cell>0;

                elseif(size(param.resolution.objects,1)>1)
                    %D. Sullivan 3/15/13, temporarily changing from max to
                    %min so that images will be at highest resolution
                    %possible, but artifacts from upsampling may be present
                    %
                    if ~isfield(param, 'outputres')
                        param.outputres = min(param.resolution.objects);
                    end
                end

                %DPS 10/20/14 - fill holes in the image. Use 4 connectivity to
                %ensure there are no holes in each slice. Sadly this requires a
                %for loop.

                if param.verbose
                    disp( 'Filling holes in framework' );
                end
                
                for i = 1:size(param.cell,3)
                    if ~param.instance.nucleus
                        param.nucleus(:,:,i) = imfill(param.nucleus(:,:,i),'holes');
                    end

                    if ~param.instance.cell
                        param.cell(:,:,i) = imfill(param.cell(:,:,i),'holes');
                    end
                end
            catch err
                getReport( err )
                %icaoberg 7/8/2013
                if ~strcmpi( param.synthesis, 'nucleus' )
                    warning(['CellOrganizer: No resolution specified for either cell or object class ',...
                        'assuming no resizing necessary. If this is incorrect unexpected results will occur.'])
                end
            end
        end

        %icaoberg 10/9/2012
        %save the synthetic nucleus to temporary folder
        if strcmpi( synthesis, 'nucleus' ) || ...
                strcmpi( synthesis, 'framework' ) || ...
                strcmpi( synthesis, 'all' )

            img = param.nucleus;
            %icaoberg 8/6/2012
            save( [param.temporary_results filesep 'image1.mat'], 'img' );
            imgs{length(imgs)+1} = img;
            clear img;
        end

        %if the user only requested a nucleus, stop method
        if strcmpi( synthesis, 'nucleus' )
            return
        end

        %save synthetic cell to temporary folder
        if strcmpi( synthesis, 'cell' ) || ...
                strcmpi( synthesis, 'framework' ) || ...
                strcmpi( synthesis, 'all' )
            img = param.cell;
            %icaoberg 8/6/2012
            save( [param.temporary_results filesep 'image2.mat'], 'img' );
            imgs{length(imgs)+1} = img;
            clear img;
        end

        %if the user requested a cell membrane or a framework, stop
        %method
        if strcmpi( synthesis, 'cell' ) || strcmpi( synthesis, 'framework' )
            return
        end

        if param.verbose
            disp('Computing Euclidean distance transform to cell and nuclear edges');
        end
        param.celldist = bwdist(bwperim(param.cell));
        param.nucleardist = bwdist(bwperim(param.nucleus));
        param.nucleardist(param.nucleus==1) = -param.nucleardist(param.nucleus==1);

        if exist( 'centrosome', 'var' ) && ~isempty( centrosome )
            fprintf( 1, '%s\n', 'Generating centrosome.' );
            [param.centrosome,outres] = model2instance( models{centrosome{1}}.proteinModel, param );
            param.resolution = ml_initparam(param.resolution,struct('centrosome',outres));
        end

        if param.verbose
            disp('Saving temporary cell shape model instance');
        end

        clear blank;
        clear box;
        clear instance;
        clear nuclei;
        % xruan 09/10/2015
        % first check the resolution setting and get full final resolution
        % in case it is a scalar to allow model2instance to use

        if( isfield( param, 'resolution' ) && isfield(param, 'outputres'))
            %since the cell is built conditionally on the nucleus, they should
            %ALWAYS have the same resolution
            initres = param.resolution.cell;%e.g. [0.05,0.05,0.35]

            %if outputres is a single value, we assume it's a scalar on the
            %image resolution
            finalres = param.outputres;

            if length(finalres) == 1
                finalres = repmat(1/finalres, [1,3]).*initres;
            end
            param.outputres = finalres;
        end
        for j=1:1:length(models)
            param.currentmodelnum = j;
            if isfield(models{j}, 'proteinModel')
                if isfield( models{j}.proteinModel, 'class' ) && ...
                        ~strcmpi( models{j}.proteinModel.class, 'centrosome' )
                    if param.verbose
                        tic
                        disp( ['Generating protein pattern model instance ' num2str(j) ' from ' models{j}.proteinModel.class ' pattern'] );
                    end

                    [img,outres] = model2instance( models{j}.proteinModel, param );
                    %img = ml_bcimg(double(img),[],[]);

                    if param.verbose
                        disp( ['Saving temporary protein model instance ' num2str(j+2) ] );
                    end
                    %icaoberg 8/6/2012
                    save( [param.temporary_results filesep 'image' num2str(j+2) '.mat'], 'img' );
                    imgs{length(imgs)+1} = img;
                    clear img;
                end
                %D. Sullivan 7/18/13 adding object movement model
                if param.randomwalk == true
                    objectMovementModel(models{j},param)
                end
            end
        end

        %D.Sullivan 1/21/13 readjust cell and nucleus to be same
        %resolution as final object models and resave temporary files

        %This code should only be relevant if there is more than one object
        %model specified. Otherwise the param.resolution.cell should have
        %been set to the same size as the objects from line 618-D Sullivan
        %10/22/14
        if( isfield( param, 'resolution' ) && isfield(param, 'outputres'))
            %since the cell is built conditionally on the nucleus, they should
            %ALWAYS have the same resolution
            % initres = param.resolution.cell;%e.g. [0.05,0.05,0.35]

            %if outputres is a single value, we assume it's a scalar on the
            %image resolution
            % finalres = param.outputres;

            % if length(finalres) == 1
            %     finalres = repmat(1/finalres, [1,3]).*initres;
            % end

            %nucleus
            origsize_x = size(param.nucleus,2);
            origsize_y = size(param.nucleus,1);
            origsize_z = size(param.nucleus,3);
            finalsize_x = ceil(initres(1)./finalres(1).*origsize_x);
            finalsize_y = ceil(initres(2)./finalres(2).*origsize_y);
            finalsize_z = ceil(initres(3)./finalres(3).*origsize_z);
            param.nucleus = imresize(param.nucleus,[finalsize_y finalsize_x],'bilinear');

            %need to resize the z
            param.nucleus = tp_stretch3d(param.nucleus,finalsize_z);
            param.nucleus = param.nucleus>0;
            imgs{1} = param.nucleus;

            if isstruct(param.nucmesh)
                param.nucmesh.vertices(:, 1) = param.nucmesh.vertices(:, 1) * origsize_x / finalsize_x;
                param.nucmesh.vertices(:, 2) = param.nucmesh.vertices(:, 2) * origsize_y / finalsize_y;
                param.nucmesh.vertices(:, 3) = param.nucmesh.vertices(:, 3) * origsize_z / finalsize_z;
            end

            %resave the synthetic nucleus to temporary folder
            if strcmpi( synthesis, 'nucleus' ) || ...
                    strcmpi( synthesis, 'framework' ) || ...
                    strcmpi( synthesis, 'all' )
                img = param.nucleus;
                %icaoberg 8/6/2012
                save( '-v7.3', [param.temporary_results filesep 'image1.mat'], 'img' );
                clear img;
            end

            %cell
            %Note: the recalculation of final size should be unnecessary since the
            %cell and nucleus images should always be the same size, but since the
            %arithmatic is trivially fast I recalculate to deal with potential
            %weird situations in the future(e.g. if we need the nuclear image to be
            %a smaller object that we add in to the cell image for space) DPS
            origsize_x = size(param.cell,2);
            origsize_y = size(param.cell,1);
            origsize_z = size(param.cell,3);
            finalsize_x = ceil(initres(1)./finalres(1).*origsize_x);
            finalsize_y = ceil(initres(2)./finalres(2).*origsize_y);
            finalsize_z = ceil(initres(3)./finalres(3).*origsize_z);
            param.cell = imresize(param.cell,[finalsize_y finalsize_x],'bilinear');

            %need to resize the z
            param.cell = tp_stretch3d(param.cell,finalsize_z);
            param.cell = param.cell;
            param.cell = param.cell>0;
            imgs{2} = param.cell;

            if isstruct(param.cellmesh)
                param.cellmesh.vertices(:, 1) = param.cellmesh.vertices(:, 1) * origsize_x / finalsize_x;
                param.cellmesh.vertices(:, 2) = param.cellmesh.vertices(:, 2) * origsize_y / finalsize_y;
                param.cellmesh.vertices(:, 3) = param.cellmesh.vertices(:, 3) * origsize_z / finalsize_z;
            end

            %resave synthetic cell to temporary folder
            if strcmpi( synthesis, 'cell' ) || ...
                    strcmpi( synthesis, 'framework' ) || ...
                    strcmpi( synthesis, 'all' )
                img = param.cell;
                %icaoberg 8/6/2012
                save( [param.temporary_results filesep 'image2.mat'], 'img' );
                clear img;
            end

            %resize the remaining images to the output resolution
            for j = 3:length(imgs)
                % xruan 09/10/2015

                % for some cases, the resolution of the proteins have been
                % adjusted, if so there is no need to crop the image.
                img = imgs{j};
                img = imresize(img,[finalsize_y finalsize_x],'bilinear');
                img = tp_stretch3d(img, finalsize_z);
                save( '-v7.3', [ param.temporary_results filesep 'image' num2str(j) '.mat'], 'img' );
                imgs{j} = img;
            end

            %icaoberg
            %super hacky must fix in next release
            microtubule_temp_folder = [ param.temporary_results filesep 'microtubules' ];
            if exist( microtubule_temp_folder )
                try
                    number_of_files = length( dir( [microtubule_temp_folder filesep '*.mat'] ) );
                catch
                    number_of_files = 0;
                end

                filename = [ microtubule_temp_folder filesep ...
                    'microtubule' num2str(number_of_files) '.mat' ];
                if exist( filename )
                    scaled_microtubules = {};
                    disp( ['Loading file ' filename] );
                    microtubules = load( filename );

                    %clear temp filename
                    for k=1:1:length(microtubules.imXYZ)
                        G = zeros( size(microtubules.img) );

                        t=1;
                        for t = 2:size(microtubules.imXYZ{k},2)
                            if prod([microtubules.imXYZ{k}(1,t), ...
                                    microtubules.imXYZ{k}(2,t), ...
                                    microtubules.imXYZ{k}(3,t)])~=0
                                G( microtubules.imXYZ{k}(1,t), ...
                                    microtubules.imXYZ{k}(2,t), ...
                                    microtubules.imXYZ{k}(3,t)) = ...
                                    G(microtubules.imXYZ{k}(1,t), ...
                                    microtubules.imXYZ{k}(2,t), ...
                                    microtubules.imXYZ{k}(3,t)) + 1;
                            end
                        end

                        low_high = [0 max(G(:))];
                        G = ml_bcimg(G,low_high,[0 255]);
                        G = imresize(G,[finalsize_y finalsize_x],'bilinear');
                        G = tp_stretch3d(G, finalsize_z);

                        [x,y,z]=ind2sub(size(G), find(G~=0));
                        scaled_microtubules{ ...
                            length(scaled_microtubules)+1} = ...
                            [x,y,z];
                        clear x y z G low_high
                    end

                    microtubules.scaled_microtubules = scaled_microtubules;
                    disp(['Updating file ' filename]);
                    save( filename, 'microtubules' );
                end
            end
        end
    otherwise
        error( 'CellOrganizer: Unknown dimensionality or dimensionality argument missing.' );
end
end%model2img
