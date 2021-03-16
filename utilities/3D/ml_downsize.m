function new = ml_downsize(old, ratio, dsmeth)
% FUNCTION NEW = ML_DOWNSIZE(OLD, RATIO, DSMETH)
% Downsize a 3D image using the downsize ratio
% Important note: this version just works for int ratio value.  Check
%     ml_3dimresize for real ratios
% old: the original image
% ratio: a length 3 vector, ratio(1) is the downsample ratio on x, ratio(2)
%        for y and ratio(3) for z.  All 3 dimensions will be halfed if a
%        [2 2 2] is used as ratio
% dsmeth: sum ('summation',blank) or average ('average'), jnewberg 11/24/05
%
% Xiang Chen
% Aug 15, 2002
%
% Aug 26, 2013 G. Johnson  Allow for insertion of 2D objects.
% Nov  4, 2013 G. Johnson  Allows for more accurate downsampling.


% Copyright (C) 2006 Murphy Lab
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

%dsmeth options are
%  'nearest' - nearest neighbor interpolation
%  'linear'  - linear interpolation
%  'spline'  - spline interpolation
%  'cubic'   - cubic interpolation as long as the data is uniformly
%              spaced, otherwise the same as 'spline'

if ~exist( 'dsmeth','var')
    dsmeth = 'linear';
end

imsize = size(old);
if length(imsize) == 2
    imsize = [imsize 1];
end

nodim = find(imsize == 1);
isdim = find(imsize ~= 1);


old = squeeze(permute(old, [isdim, nodim]));
ratio = ratio(isdim);
imsize = imsize(isdim);

% if length(ratio) == 1
%     ratio = repmat(ratio, [1,3]);
% elseif length(ratio) == 2;
%     ratio = [ratio 1];
% end

c = class(old);



%conver to double to linearly interpolate
old = double(old);

dim = round(imsize./ratio); %% desired output dimensions rounded to nearest pixel

argin = cell(size(ratio));
for i = 1:length(ratio)
    argin{i} = linspace(1, imsize(i), dim(i));
end

[interps{1:length(argin)}] = ndgrid(argin{:});
    
      
if strcmpi(dsmeth, 'nearest')
    for i = 1:length(interps)
        interps{i} = round(interps{i});
    end
end

new = interpn(old, interps{:}, dsmeth, 0);

new = permute(new, [isdim, nodim]);

new = eval([ c '(new)']);
end