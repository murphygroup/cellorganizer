function centralmoment = ml_imgcentmoments(image, x, y)
% ML_IMGCENTMOMENTS(IMAGE, X, Y) calculates the central moment MUxy for IMAGE
% ML_IMGCENTMOMENTS(IMAGE, X, Y), 
%    where IMAGE is the image to be processed and X and Y define
%    the order of the moment to be calculated. The coordinate system 
%    is centered using the center of fluorescence (m10/m00, m01/m00).
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

% 19 Aug 98 - M.V. Boland

% $Id: ml_imgcentmoments.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

if nargin ~= 3
    error('Please supply all three arguments (IMAGE, X, Y)') ;
end

if (isempty(image))
    error('IMAGE is empty.') 
end

% 
% Convert image to double precision if necessary
%
if (~isa(image,'double'))
    image = double(image) ;
end

m00 = ml_imgmoments(image, 0, 0) ;
m10 = ml_imgmoments(image, 1, 0) ;
m01 = ml_imgmoments(image, 0, 1) ;

cofx = m10/m00 ;
cofy = m01/m00 ;

%
% Generate a matrix with the x coordinates of each pixel.
%  If the order of the moment in x is 0, then generate
%  a matrix of ones.  Subtract the center of fluorescence
%  from the x coordinates.
%
if x==0
    xcoords = ones(size(image)) ;
else
    xcoords = (ones(size(image,1),1) * (([1:size(image,2)] - cofx) .^ x)) ;
end

%
% Generate a matrix with the y coordinates of each pixel.
%  If the order of the moment in y is 0, then generate
%  a matrix of ones Subtract the center of fluorescence 
%  from the y coordinates.
%
if y==0
    ycoords = ones(size(image)) ;
else
    ycoords = ((([1:size(image,1)]' - cofy) .^ y) * ones(1,size(image,2)));
end

%
% Multiply the x and y coordinate values together
%
xycoords = xcoords .* ycoords ;

%
% The central moment is the double sum of the xyf(x,y)
%
centralmoment = sum(sum(xycoords .* image)) ;
