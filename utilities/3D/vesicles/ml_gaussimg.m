function img = ml_gaussimg( sigma, param )
%ML_GAUSSIMG Synthesize an image from a 3D Gaussian distribution.
%   IMG = ML_GAUSSIMG(SIGMA) returns an image that has intensities with
%   a 3D Gaussian distribution. SIGMA is the 2x2 covariance matrix of the 
%   Gaussian distribution.

% 26-Jan-2006 T. Zhao
%
% Copyright (C) 2007-2012 Murphy Lab
% Carnegie Mellon University
%
% August 4, 2012 D. Sullivan added param structure to the method, param.imagesize
%       contains the size of the cell image which the objects will be added.
%       This is used to make sure we don't waste memory making object
%       images that can't fit in the final image. 
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

%gj 2/15/2013 Set param as optional
%gj 2/15/2013 Because of numerical error, if we have a covariance matrix
%with very small values in the Z direction, the gaussian with sum to >>> 1.
%If the image is 1 pixel deep, we enforces sum(img(:)) == 1

if nargin < 1
    error('Exactly 1 arguments are required')
end

ndims = size(sigma,1);
zeropos = zeros(1,ndims);


%8/15/13 D. Sullivan Added check for param.objstd which indicates the
%number of standard deviations at which to synthesize objects.
if ~exist('param', 'var') || ~isfield(param,'objstd')
%     warning('No standard deviation specified for object synthesis, defaulting to 6.');
    param.objstd = 6;
end

 [a,e] = eig(sigma);
 imsize = ceil((sqrt(diag(sigma)).*param.objstd).*2)';

% imsize = ceil(param.objstd * sqrt(diag(sigma)))';

%devins 8/4/2012 check if too large for image
if exist('param', 'var') & isfield(param, 'imagesize')
    imsize = min([imsize;param.imagesize]);
end

x = ml_imcoords(imsize,1,-round(imsize/2))';
img = reshape(mvnpdf(x,zeropos,sigma),imsize);

if sum(img(:)) ~= 1
    img = img / sum(img(:));
end

if exist('param', 'var') & isfield(param, 'rendAtStd')
    zeropos = zeros(size(sigma,1),1);
    pos = zeros(size(sigma,1),1);
    
    [a,e] = eig(sigma);
    pos(1) = param.rendAtStd*sqrt(e(1,1));
    img(img <= mvnpdf(pos, zeropos, e)) = 0;
end
    
end