function moment = ml_imgmoments(image, x, y)
% ML_IMGMOMENTS(IMAGE, X, Y) calculates the moment MXY for IMAGE
% ML_IMGMOMENTS(IMAGE, X, Y), 
%    where IMAGE is the image to be processed and X and Y define
%    the order of the moment to be calculated. For example, 
%    ML_IMGMOMENTS(IMAGE,0,1) calculates the first order moment 
%    in the y-direction, and 
%    ML_IMGMOMENTS(IMAGE,0,1)/ML_IMGMOMENTS(IMAGE,0,0) is the 
%    'center of mass (fluorescence)' in the y-direction
%

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

% 10 Aug 98 - M.V. Boland

% $Id: ml_imgmoments.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

if nargin ~= 3
    error('Please supply all three arguments (IMAGE, X, Y)') ;
end

%
% Check for a valid image and convert to double precision
%   if necessary.
%
if (isempty(image))
    error('IMAGE is empty.') 
elseif (~isa(image,'double'))
    image = double(image) ;
end

%
% Generate a matrix with the x coordinates of each pixel.
%  If the order of the moment in x is 0, then generate
%  a matrix of ones
%
if x==0
    if y==0
        xcoords = ones(size(image)) ;
    end
else 
    xcoords = (ones(size(image,1),1) * ([1:size(image,2)] .^ x)) ;
end

%
% Generate a matrix with the y coordinates of each pixel.
%  If the order of the moment in y is 0, then generate
%  a matrix of ones
%
if y~=0
    %	ycoords = ones(size(image)) ;
    ycoords = (([1:size(image,1)]' .^ y) * ones(1,size(image,2))) ;
end

%
% Multiply the x and y coordinate values together
%
if y==0
    xycoords = xcoords ;
elseif x==0
    xycoords = ycoords ;
else
    xycoords = xcoords .* ycoords ;
end

%
% The moment is the double sum of the xyf(x,y)
%
moment = sum(sum(xycoords .* image)) ;
