function [isub] = ml_imgbgsub(image, method)
%  ML_IMGBGSUB - Subtract the background from an image
%  [ISUB] = ML_IMGBGSUB(IMAGE, METHOD)
%
%    Outputs:
%     ISUB - IMAGE with the background subtracted
%
%    Inputs:
%     IMAGE - image from which background should be subtracted
%     METHOD - method to use for identifying background
%              'common' - use the most common pixel value
%              'lowcommon' - the most common pixel between 0 and the mean
%                       of the image
%              'quantile' - the most common pixel between 0 and 0.75 
%                       quantile of the image.

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

%    M. Boland - 01 Mar 1999
%    05 Nov 2005 T. Zhao add 'lowcommon' method
%    07 Nov 2005 T. Zhao add 'quantile' method

% $Id: ml_imgbgsub.m,v 1.6 2006/06/27 13:33:47 tingz Exp $

if (isempty(image)),
    error('IMAGE is empty') ;
end

methods = {'common','lowcommon','quantile'} ;
if sum(strcmp(method, methods)) == 0
    error('Undefined method for determining the most common pixel');
end


% if (strcmp(methods,'common')),
% 'lowcommon', 'quantile' methods added by T. Zhao  
switch method
case 'common'
    %
    % Find the most common pixel value
    common = ml_imgcommonpixel(image) ;
case 'lowcommon'
    common = ml_imgcommonpixel(image,mean(image(:))-1) ;
case 'quantile'
    common = ml_imgcommonpixel(image,prctile(image(:),75)-1) ;
end 

%
% Check for a common pixel value of 0
if (common == 0),
    isub = image ;
    return
end

%t+ Feb. 23, 2006
image = double(image);
%t++

%
% Subtract the most common pixel value from each pixel in the image
isub = image - common ;

%
% Set any pixel values < 0 to 0
isub = (isub>0).*isub ;

%
% endif

