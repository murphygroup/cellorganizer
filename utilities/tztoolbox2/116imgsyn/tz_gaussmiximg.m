function img = tz_gaussmiximg(mix)
%TZ_GAUSSMIXIMG Obsolete. See ML_GAUSSMIXIMG.
%   IMG = TZ_GAUSSMIXIMG(MIX) returns an image that has intensities with
%   a 2D Gaussian mixture distribution specified by the structure MIX. See
%   GMM for details about the structure.
%   
%   See also GMM GMMPROB

%   10-Jul-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu

error(tz_genmsg('of','tz_gaussmiximg','ml_gaussmiximg'));

if nargin < 1
    error('Exactly 1 argument is required');
end

boxes = [];
for i=1:mix.ncentres
    boxSize(i,:) = [6*sqrt(mix.covars(1,1,i)),6*sqrt(mix.covars(2,2,i))];
    boxes = [boxes;mix.centres(i,:)-boxSize(i,:)/2; ...
        mix.centres(i,:)+boxSize(i,:)/2];
end
    
topLeft = min(boxes,[],1);
bottormRight = max(boxes,[],1);
imageSize = round(bottormRight-topLeft)+[1 1];

offset = topLeft;

x = tz_imcoords(imageSize,1,offset-1)';
% mix.centres = ml_addrow(mix.centres,offset);

img = reshape(gmmprob(mix,x),imageSize(1),imageSize(2));