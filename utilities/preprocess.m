function [segimg, top_slice, bot_slice, sumslices, zflip ] = preprocess( img, varargin )
%PREPROCESS segments and preprocesses a given image.
%
% Backwards compatibility mode
% Inputs
% ------
% img = image array to be processed
% imgFile = mat file containing image to be processed
% psfPath = path to psf to use
% downsample = [1xn] or single number defining the amount to downsample the image before preprocessing.
% adherent = boolean flag defining if the cell is adherent(forces largest slice to bottom)
% top_thresh = the fraction of fluorescence you wish to consider as the top slice
% bot_thresh = the fraction of fluorescence you wish to consider as the bottom slice
% display = boolean flag of whether to display progress on screen
% cellmask = binary valued image array masking the inside of the cell
%
% Outputs
% -------
% segimg = segmented image
% top_slice = integer value of the top slice where signal lies
% bot_slice = integer value of the top slice where signal lies
% sumslices = the cumsum of each slice
% zflip = boolean flag indicating if preprocess thinks the cell is upsidedown

% Devin Sullivan
%
% Copyright (C) 2006-2016 Murphy Lab
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

% Created by Devin Sullivan 3/21/13 moved from being a subfunction of
% seg_cell_and_nucleus.m
%
% G. Johnson 8/22/13 - support for 2D masks, even if img is 3D
%
% G. Johnson 8/26/13 - Re-mask the output of region_seg incase the snake moves
%                     out of bounds
%
% G. Johnson 1/15/14 - Remove psf as per Ivan's request
%
% D. Sullivan 4/24/14 Removed if statment that would not work if bot_slice and bot_seg were
%                     equal
%
%icaoberg 26/1/2014
%backward compatibility mode means that it will run as before but new call
%will simplify calling this method directly
%
% R.F. Murphy 3/12/2017 - add active contour parameters specific to preprocessing
%                         to param struct (keeping same defaults) and 
%                         cleanup arg processing 

if nargin == 9 | nargin == 10
    disp( 'Running method in backward compatibility mode' )
    imgFile = varargin{1};
    psfPath = varargin{2};
    downsample = varargin{3};
    adherent = varargin{4};
    top_thresh = varargin{5};
    bot_thresh = varargin{6};
    display = varargin{7};
    cellmask = varargin{8};
    param = [];
    if nargin == 10
        param = varargin{9};
    end
elseif nargin == 2
    if ~isfield( param, 'image_file' )
        imgFile = '';
    else
        imgFile = param.image_file;
    end
    
    if ~isfield( param, 'psf_path' );
        psfPath = '';
    else
        psfPath = varargin{2};
    end
    
    if ~isfield( param, 'downsample' )
        downsample = [ 5, 5, 1 ];
    else
        downsample = param.downsample;
    end
    
    if ~isfield( param, 'adherent' )
        adherent = true;
    else
        adherent = param.adherent;
    end
    
    if ~isfield( param, 'top_threshold' )
        top_thresh = 0.98;
    else
        top_thresh = param.top_threshold;
    end
    
    if ~isfield( param, 'bottom_threshold' )
        bot_thresh = 0.02;
    else
        bot_thresh = param.bottom_threshold;
    end
    
    if ~isfield( param, 'cellmask' )
        cellmask = [];
    else
        cellmask = varargin{8};
    end
    
    param = [];
else
    warning( 'Wrong number of input arguments. Exiting method.' );
    segimg = []; top_slice = []; bot_slice = []; sumslices = [];  zflip = [];
    return
end

param = ml_initparam( param, struct( ...
    'display', false, ...
    'debug', false, ...
    'verbose', false, ...
    'preprocess', [] ) );
% note that there is also param.stiffness which controls nuclear hole finding
param.preprocess = ml_initparam( param.preprocess, struct( ...
    'stiffness', 0.7, ...
    'maxiter', 3000, ...
    'quittol', 0.00001 ) );

if ~exist('cellmask', 'var') || isempty( cellmask )
    cellmask = ones(size(img));
end

imsize = size(img);
masksize = size(cellmask);

%icaoberg
%if the psfpath is empty then do not attempt to load psf
% Load PSFs
if ~isempty( psfPath )
    infPSF = imfinfo( psfPath);
    
    psf = zeros( infPSF(1).Height, infPSF(1).Width, length(infPSF));
    
    for I=1:length(infPSF)
        psf(:,:,I)=imread(psfPath,I);
    end
    
    psf = psf.^2; % Approximate confocal PSF
else
    psf = [];
end

%D. Sullivan 3/21/13 initialize the z flip to 0
zflip = 0;

%icaoberg 15/1/2014
%downsample PSFs and images as specified by user
%if the psf is not specified then the psf should be empty
if exist( imgFile, 'file' )
    load(imgFile)
else
    if ~isempty( psf )
        if param.verbose
            disp( 'Downsampling PSF and image');
            tic;
        end
    else
        if param.verbose
            disp( 'Downsampling image');
            tic;
        end
    end
    
    if param.verbose
        disp( [ 'Downsampling using vector [x,y,z]: [' ...
            num2str(downsample) ']' ] );
        disp( [ 'Image size [' num2str(size(img) ) ']' ] );
        try
            disp( [ 'Image resolution [' ...
            num2str(param.model.resolution) ']' ] );
        catch
            disp( ['Image resolution not found.'] );
        end
    end
    
    if ~isempty( psf )
        psf = ml_downsize(psf,downsample,'linear');
    end
    
    resizeimg = ml_downsize( img, downsample, 'linear'); 
    
    if param.verbose
    toc
    end
    
    %     if length(size(cellmask)) == 3 && masksize(3) ~= imsize(3)
    %         [~,ind] = max(squeeze(sum(sum(cellmask,1),2)));
    %         cellmask = cellmask(:,:,ind);
    %     end
    
    if length(size(cellmask)) == 3
        %if the mask isnt the same size as the currently downsampled image
        if ~all(size(cellmask) == size(resizeimg))
            resizecellmask = ml_downsize(cellmask,downsample) > 0;
        else
            resizecellmask = cellmask;
        end
        
    else
        %D. Sullivan 5/5/14 - added check to see if downsampling is
        %necessary
        %if the mask isnt the same size as the currently downsampled image
        if (size(cellmask,1) ~= size(resizeimg,1)) || (size(cellmask,2)~=size(resizeimg,2))
            resizecellmask = ml_downsize(cellmask, downsample(1:2)) > 0;
        else
            resizecellmask = cellmask;
        end
        resizecellmask = repmat(resizecellmask, [1,1, size(resizeimg,3)]);
    end
    
    clear dnaim3
    
    %icaoberg 15/1/2014
    %only attempt to deconvolve image if psf is not empty
    if ~isempty( psf )
        if param.verbose
            disp('Deconvolving image');
            tic;
        end
        
        [image,psf2] = deconvblind(resizeimg,psf);
        
        if param.verbose
            toc
        end
    else
        image = resizeimg;
    end
end

if param.verbose
    disp('Segmenting image');
    tic;
end

sumslices = cumsum(squeeze(sum(sum(image))));
sumslices = sumslices/sumslices(end);
bot_slice=find(sumslices> bot_thresh); %0.05, 0.02
bot_slice=bot_slice(1);
top_slice=find(sumslices< top_thresh); %0.95, 0.98
%D. Sullivan 5/5/14 - if there are no slices less than the top_thresh, set
%top_slice to maximum value in image
if isempty(top_slice)
    top_slice = size(image,3);
else
    top_slice=top_slice(end);
end
%D.Sullivan 6/6/13, these two lines literally do nothing
% bot_slice=max(bot_slice,bot_slice);
% top_slice=min(top_slice,top_slice);

if param.verbose
    disp( ['Total number of slices: ' num2str(length(sumslices))] );
    disp( ['Bottom slice index: ' num2str(bot_slice)] );
    disp( ['Top slice index: ' num2str(top_slice)] );
end

mask = image > ml_rcthreshold(image(:));
for i = 1:size(mask,3)
    mask(:,:,i) = bwfill(mask(:,:,i), 'holes');
end
mask = repmat(sum(mask,3), [1,1,size(mask,3)]) > 0;

% bounds = ones(size(mask));
% bounds(3:end-3,3:end-3,:) = 0;
% mask(bounds >0) = 0;

if param.verbose
    disp( 'Cropping 3D image using mask' ); 
    tic
end

mask = tz_maskimg_3d( mask, resizecellmask );
    
if param.verbose
    toc
end

if param.verbose
    disp( 'Region Based Active Contour Segmentation' ); 
    tic
end

alpha = param.preprocess.stiffness;
display = param.display;
maximum_iterations = param.preprocess.maxiter;
quit_tolerance = param.preprocess.quittol;

%icaoberg 26/1/2014
%pass in parameter structure
segimg = region_seg(image, mask, ...
    maximum_iterations, alpha, display, quit_tolerance, param );
if param.verbose
    toc
end

clear directory
clear directory2
clear output_filename

segimg = ml_findmainobj(segimg);
segimg = tz_maskimg_3d(segimg, resizecellmask);

inds = find(sum(sum(segimg,1),2));
bot_seg = inds(1);
top_seg = inds(end);

bot_slice=max([bot_slice,bot_seg]);
top_slice=min([top_slice,top_seg]);

%devins 24/4/2014
%removed if statment that would not work if bot_slice and bot_seg were
%equal
if adherent %false, <param>
    %D. Sullivan  3/21/13 check if image is right side up.
    %since we expect the cell to generally get larger at the bottom for
    %adherent cell lines, if the slope of the areas is positive we should
    %flip the order
    
    %get total cell area per slice
    areas = squeeze(sum(sum(segimg)));
    %eliminate ones that are not in the cell
    cellareas = areas(areas~=0);
    %P = coefficients P(1)*X+P(2)
    P = polyfit([1:length(cellareas)],cellareas',1);
    if P(1)>0
        %set zflip flag
        zflip = 1;
        segimg = flipdim(segimg,3);
    end
    
    [bot_slice,top_slice]=fixadherent(segimg,bot_slice,top_slice);
end
clear image