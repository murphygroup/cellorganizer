function out_img = im2projection( img, param )
% IM2GPROJECTION creates a sum or mean projection of the input image
%
% List Of Input Arguments   Descriptions
% -----------------------   ------------
% img                       3D binary or realvalued image.
% param                     struct with a 'method' field that can be set
%                           to 'mean' if a mean value projection is desired
%
% List Of Outputs     Descriptions
% ---------------     ------------
% out_img             a 2D image that contains a projection
%                     in each dimension of the original image

% Author: Yue Yu and Ivan Cao-Berg
%
% Copyright (C) 2012-2016 Murphy Lab
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

if nargin == 1
    param.method = 'mean';
end

out_img = [];
if isempty( img )
    warning( 'Input argument image is empty' );
    return
end

try
    [m,n,q] = size(img);
    out_img = zeros((m+q),(q+n));
    
    if strcmpi( param.method, 'sum' )
        sumq = sum(img,3);
        summ = squeeze(sum(img,1));
        sumn = squeeze(sum(img,2));
        out_img(1:m,q+1:end) = sumq;
        out_img(1:m,1:q) = sumn;
        out_img(m+1:end,q+1:end) = flipud(summ');
    elseif strcmpi( param.method, 'mean' )
        sumq = mean(img,3);
        summ = squeeze(mean(img,1));
        sumn = squeeze(mean(img,2));
        out_img(1:m,q+1:end) = sumq;
        out_img(1:m,1:q) = sumn;
        out_img(m+1:end,q+1:end) = flipud(summ');
    else
        out_img = [];
        warning('Unknown method');
    end
catch err
    out_img = [];
    warning('Unable to make projections');
    getReport( err )
    return
end
