function [procimage,imagemask,resimg] = ...
    ml_preprocess( image, cropimage, way, bgsub,threshmeth,highpass)

% [PROCIMAGE,IMAGEMASK,RESIMG] = ML_PREPROCESS( IMAGE, CROPIMAGE, WAY, BGSUB)
%
% Performs background subtraction, cropping with CROPIMAGE, and
% thresholding exactly in the said order. This is the standard
% Murphy Lab way of doing it. Note that in Michael's feature code
% the order of cropping and thesholding is reversed.
% If you have no crop image, give CROPIMAGE as [].
% PROCIMAGE is the processed image ready
% to have its features calculated by ml_features()
% IMAGEMASK is the binary mask of procimage with noise removed
% by majority filtering.
% RESIMG is the residual image left over after thresholding, i.e.
% the below-threshold image (but it is background subtracted).
% WAY determines whether preprocessing should be done Michael
% Boland way (threshold before crop), or Murphy Lab way (crop
% before threshold).
% 'ml' = Murphy Lab way
% 'mb' = Michael Boland way
% BGSUB should be 'nobgsub' if you do not want background
% subtraction, otherwise omit this argument. For more details, see
% ml_imgbgsub
%
%Edited: 
%6/20/13 D. Sullivan - added highpass image filtering
  
% Copyright (C) 2006  Murphy Lab
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

%fprintf(1,'ml_preprocess being called...\n');

if( nargin < 4)
    bgsub = 'yesbgsub';
end

if nargin<5
    threshmeth = 'nih';
end

if nargin==6 && highpass==1
    %D. Sullivan 6/20/13 - this should be resolution and/or image size
    %dependent in the future
    filtersize = 11;
    I_filt = fspecial('disk',filtersize);
    blurred_image = imfilter(image,I_filt,'same',mean(image(:)));
    image = image - blurred_image;
    %recontrast stretch it to prevent negative numbers 
    image = uint8((image-min(image(:)))./max(image(:)).*255);
    
end

switch bgsub
case 'nobgsub'
otherwise
    if strcmp(bgsub,'yesbgsub')
        bgsub = 'common';
    end
    %%% subtract background
    image = ml_imgbgsub( image, bgsub );
end

if( isempty( cropimage))
    threshcropimage = [];
else
    switch way
    case 'ml'
        %%% Crop (or mask rather)
        image( find( cropimage==0)) = 0;
        % Make sure no cropping is done later
        threshcropimage = [];
    case 'mb'
        % Make sure cropping is done later
        threshcropimage = cropimage;
    otherwise
        error('Unrecognized option for the WAY argument');
    end
end

if any(image(:)<0)
    error('The intensity of the image should not be negative.');
end

%%% Theshold
if (max(image(:))>0)
    if strcmp(threshmeth,'none')
        imagemask = image>0;
        procimage = image;
        resimg = zeros(size(image));
    else
	scale = 2;
        imagemask = ml_threshcrop(image, threshcropimage,threshmeth, scale);

        if (max(imagemask(:))>0)
%             procimage = roifilt2(0, image, ~imagemask);
            procimage = image;
            procimage(find(imagemask==0)) = 0;
            resimg = image;
            resimg(find(imagemask))=0;
        else
            warning('No protein fluorescence left after preprocessing');
            resimg = [];
            procimage = [];
        end
    end
else
    warning('No protein fluorescence left after preprocessing');
    imagemask = [];
    resimg = [];
    procimage = [];
end

