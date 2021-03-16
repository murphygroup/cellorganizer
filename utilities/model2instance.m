function [object, outres] = model2instance( model, param )
% MODEL2INSTANCE Generate a single instance from a valid SLML model.
% Model is a valid protein pattern model. Depending on the model type,
% the structure of param changes.

% Author: Ivan E. Cao-Berg (icaoberg@scs.cmu.edu)
%
% Copyright (C) 2011-2016 Murphy Lab
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

% March 8, 2012 Added location option to param structure and set a default value to 'all'
%
% March 21, 2012 Changed protein.location parameter to protein.cytonuclearflag
%
% October 17, 2012 D. Sullivan added param.cellresolution and param.objresolution
% and resizing code to make everything the resolution of the object model
% being synthesized.
%
% November 14, 2012 D. Sullivan fixed bugs in resizing cell/nuclear images
% to match object size. Needed to reset param.cell and param.nucleus and
% resize the param.nucleardist and param.celldist images as well
%
% November 15, 2012 I. Cao-Berg Changed warning to regular display
%
% January 21, 2013 D. Sullivan updated resolution framework s.t. user may
% now specify multiple object model resolutions and the output will be in
% the form of the lowest resolution (highest numbers since resolutions =
% microns/pixel)
%
% January 21, 2013 D. Sullivan fixed bug where file was returning if no PSF
% was used. This was necessary to have the file run to the end and resize
% to the output resolution. Moved PSF application to after image resize.
%
% January 21, 2013 D. Sullivan added resolution field to psf and resized
% based on output resolution. Prints message to user informing of resize.
%
% January 22, 2013 I. Cao-Berg modified if statement so that it checks
% whether field exists before querying it
%
% February 4, 2013 D. Sullivan added outres argument to arguments returned.
%
% February 9, 2013 D. Sullivan added imXYZ which is a cell array of each MT
% position in addition to the outputs already returned by
% model2microtubules
%
% February 20, 2013 D. Sullivan fixed bug which was forcing 'sampled'
% version of the code to be binary.
%
% February 27, 2013 D. Sullivan fixed up/downsampling for all patterns to
% occur in this function so that it is done consistently and functions
% downsampling things only upsample again after passing back a smaller data
% matrix (should speed things up a bit)
%
% September 25, 2013 D. Sullivan added corrections for rounding errors that
% allowed objects that were resized to end up outside the cell.
%
% May 25, 2016 I. Cao-Berg Fixed references to model class and type
%
% November 13, 2016 I. Cao-Berg Updated switch to match the model classes
% and types described in the documentation
%
% January 24, 2017 Xiongtao Ruan Added reference to T cell model
% 
% February 5, 2017 I. Cao-Berg Modified function so that microtubules patterns
%                  remain grayscale after adjusting resolution
% 
% May 21, 2019 Xiongtao Ruan Added matching of resolutions of mesh


object = [];
outres=[];

if nargin > 2
    error('CellOrganizer: Wrong number of input arguments.');
end

if length(model) > 1
    if isstruct(model)
        model=model(1);
    else
        model = model{1};
    end
end

try
    dimensionality = model.dimensionality;
catch
    warning( 'CellOrganizer: Unable to set dimensionality' );
    return
end

try
    nucleus = param.nucleus;
catch
    warning('CellOrganizer: Parameter nuclear shape instance not set.')
    return
end

try
    cellMembrane = param.cell;
catch
    warning('CellOrganizer: Parameter cell shape instance not set.')
    return
end

try
    modelClass = model.class;
catch
    warning('CellOrganizer: Parameter model class not set.')
    return
end

try
    modelType = model.type;
catch
    warning('CellOrganizer: Parameter model type not set.');
    return
end

try
    microscope = param.microscope;
catch
    microscope = 'none';
end

switch lower(dimensionality)
    case '2d'
        switch lower(modelType)
            case 'vesicle'
                param = ml_initparam(param, ...
                    struct('imageSize',[1024 1024],'gentex',0,'loc','all'));
                object = ml_genprotimg( model, nucleus, cellMembrane, param );
            otherwise
                warning('CellOrganizer: Unknown or unsupported model type.');
                object = [];
                return
        end
    case '3d'
        %D. Sullivan 10/17/12
        %resample cell and nuc to proper resolution these values are expected
        %to be in um/pixel and have a value for x,y,and z
        %D. Sullivan 11/14/12
        %fixed bugs by rewriting param.cell and param.nucleus which are actually
        %what is passed to tp_genprotimage. For some reason need to continue to
        %actually resize the read in images from the param structure for the pass
        %to model2microtubules. Also needed to resize the param.nucleardist and
        %param.celldist images appropriately
        try
            cellMemOriginal = cellMembrane;
            nucleusOriginal = nucleus;
            
            %D. Sullivan 11/3/14 - select the current object resolution
            %model
            param.resolution.objects = param.resolution.objects(param.currentmodelnum,:);
        catch
            %icaoberg 11/15/2012
            %warning(['CellOrganizer: No resolution specified for either cell or object class',...
            %' assuming no resizing necessary. If this is incorrect unexpected results will occur'])
            disp(['No resolution specified for either cell or object class', ...
                ' assuming no resizing necessary. If this is incorrect unexpected results will occur'])
        end
        
        switch lower(modelType)
            case 'gmm'
                if strcmpi( model.class, 'vesicle' )
                    %tries to retrieve sampling parameters which if they exist, must already be valid
                    try
                        sampling = param.sampling.method;
                    catch
                        param.sampling.method = 0;
                    end
                    
                    try
                        N = param.sampling.N;
                    catch
                        param.sampling.N = [];
                    end
                    
                    try
                        location = param.sampling.location;
                    catch
                        param.sampling.location = 'all';
                    end
                    
                    %11/3/14 - D. Sullivan
                    %Need to add support for synthesis since it only really
                    %works in cubic voxels due to the random rotation. I'm
                    %putting an option to over
                    if isfield(param,'cubicOverride') && param.cubicOverride
                        warning( ['You have chosen to synthesize a vesicular',...
                            ' instance from non-cubic voxels. This is not recommended',...
                            ' due to random rotations imparted on the vesicles'] );
                    else
                        param.resolution.cubic = repmat(min(param.resolution.objects),1,3);
                        param.cell = AdjustResolutions(param.cell,param.resolution.cell,param.resolution.cubic);
                        param.nucleus = AdjustResolutions(param.nucleus,param.resolution.cell,param.resolution.cubic);
                        
                        % xruan 05/21/2019 change resolution for mesh
                        if isfield(param, 'nucmesh')
                            param.nucmesh.vertices = param.nucmesh.vertices .* repmat(param.resolution.cell ./ param.resolution.cubic, size(param.nucmesh.vertices, 1), 1);
                        end
                        if isfield(param, 'cellmesh')
                            param.cellmesh.vertices = param.cellmesh.vertices .* repmat(param.resolution.cell ./ param.resolution.cubic, size(param.cellmesh.vertices, 1), 1);
                        end

                        %We will also need to re-comput the distance images
                        %now - copy pasta from model2img
                        %                     fprintf( 1, '%s\n', ['Computing Euclidean distance transform to cell and nuclear edges'] );
                        param.celldist = bwdist(bwperim(param.cell));
                        param.nucleardist = bwdist(bwperim(param.nucleus));
                        param.nucleardist(param.nucleus==1) = -param.nucleardist(param.nucleus==1);
                    end
                    
                    
                    object = tp_genprotimage( model, param );
                    %D. Sullivan 1/21/13 fixed if-return so that file
                    %executes until the end
                    %D. Sullivan 1/21/13 moved psf application to after final
                    %image is created
                end
                
                if strcmpi( model.class, 'centrosome' )
                    [object,resolution] = model2centrosome(nucleus, cellMembrane, model,param);
                    param.outputres = resolution;
                end
            case 'microtubule_growth'
                try
                    centrosome = param.centrosome;
                catch
                    warning(['CellOrganizer: A centrosome instance is needed to ' ...
                        'synthesize a cytoskeletal instance. Remember to initialize param.centrosome.'])
                    object = [];
                    return;
                end
                
                %D. Sullivan 2/4/13 added resolution parameter in return
                %D. Sullivan 2/9/13 added individual MT's in return
                [object,imXYZ,resolution]  = model2microtubules( cellMembrane, nucleus, centrosome, model, param );
                %param.cubicOverride = 1;
                param.resolution.cubic = resolution;
                %D. Sullivan 1/21/13 fixed if-return so that file
                %executes until the end
                %D. Sullivan 1/21/13 moved psf application to after final
                %image is created
            
            % xruan 09/17/2016
            case 'standardized_map_half-ellipsoid'             
                [object,resolution] = model2tcell(nucleus, cellMembrane, model,param);
                param.outputres = resolution;
                
            otherwise
                warning('CellOrganizer: Unknown or unsupported model type.');
                object = [];
                return
        end
        
        %dpsull 1/21/13
        %need to resize the image to uniform final output resolution given
        %by param.outputres
        %11/3/14 - D. Sullivan
        %Need to re-adjust the image for the final output resolution.

        if (isfield(param,'cubicOverride') && param.cubicOverride) || ...
                strcmpi( modelClass, 'centrosome' ) || strcmpi( modelType, 'standardized_map_half-ellipsoid' )
        else
            %Don't need to readjust cell/nuc since they are not being
            %returned here.
            %             param.cell = AdjustResolutions(param.cell,param.resolution.cubic,param.outputres);
            %             param.nucleus = AdjustResolutions(param.nucleus,param.resolution.cubic,param.outputres);
            %Need to specify ceiling for re-adjustment becasue the previous
            %rescaling was floored
            ceilfloor = 'ceil';
            
            %icaoberg 20170204 had to change the isbinary flag for the next
            %method to false
            object = AdjustResolutions(object,param.resolution.cubic, ...
                param.outputres,false,ceilfloor);
        end
        
        %D. Sullivan 2/20/13 need to output the object as binary if
        %'disc' mode
        if strcmpi(param.sampling.method,'disc') && ...
                strcmpi(modelClass,'vesicle')
%                 strcmpi(modelType,'vesicle')
            object = object>0;
        end

        %D. Sullivan 2/4/13 add outres returning the resolution of the
        %image generated
        outres = param.outputres;
        
        %D. Sullivan 1/21/13 moved psf application to after resize.
        [psf,psfres] = micro2psf( microscope );
        %if ~isempty( psf ) && ~strcmpi(modelType,'centrosome')
        if ~isempty( psf ) && ~strcmpi(modelClass,'centrosome')
            %D. Sullivan 1/21/13
            %check if psf resolution is the same as param.outputres
            %if it doesn't resize psfres and print a message.
            if isfield(param,'outputres') && sum(psfres~=param.outputres)~=0
                disp(['Different final resolution specified than trained PSF, ',...
                    'resizing PSF to match final resolution']);
                
                outpsf_x = floor(psfres(1)./param.outputres(1).*size(psf,1));
                outpsf_y = floor(psfres(2)./param.outputres(2).*size(psf,2));
                outpsf_z = floor(psfres(3)./param.outputres(3).*size(psf,3));
                psf = imresize(psf,[outpsf_x outpsf_y],'bilinear');
                
                %need to resize the z
                psf = tp_stretch3d(psf,outpsf_z);
                
            end
            %end 1/21/13 addition
            
            object = psf_blur_hela_mean_dps( object, psf );
        end
    otherwise
        warning('CellOrganizer: Unknown or unsupported dimensionality.')
        object = [];
        return
end
