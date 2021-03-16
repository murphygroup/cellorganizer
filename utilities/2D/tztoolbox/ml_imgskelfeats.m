function [names, values, slfnames] = ml_imgskelfeats(imageproc )
% [NAMES, VALUES, SLFNAMES] = ML_IMGFEATURES(IMAGEPROC) calculates 
%    skeleton features for IMAGEPROC
%
%    where IMAGEPROC contains the pre-processed fluorescence image, 
%    Pre-processed means that the image has been cropped and had 
%    pixels of interest selected (via a threshold, for instance).
%
% Created by MV 1/18/02
% Modified MV 6/2/02: Added SLF names

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


% Initialize the variable that will contain the
% values of the features.
values = [] ;

% Find objects in the image
%
imagelabeled = bwlabel(im2bw(imageproc)) ;
obj_number = max(imagelabeled(:)) ;

% For each object ...
for (i=1:obj_number)
    % Get an image of the single object
    objmask = (imagelabeled==i) ;
    objimage = imageproc .* objmask ;
    % Compute skeleton features
    [skelfeats, names] = ml_objskelfeats( objimage);
    values = [values; skelfeats];
end

% Average the skeleton features over the whole cell
values = mean(values,1);
% SLF names
slfnames = {};
for feat_no = 1:5
    slfnames{feat_no}=['SLF7.' num2str(feat_no+79)];
end
