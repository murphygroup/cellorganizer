function [pencpixel, penflu] = ml_3dedgefeatures(protimg, protbin)

%ML_3DEDGEFEATURES calculate edge features
%   [PENCPIXEL, PENFLU] = ML_3DEDGEFEATURES(IMG, MASK) returns
%   two edge features for the 3d image IMG. MASK is the binary
%   image from IMG. Edge detection is done on MASK slice by slice.
%   PENCPIXEL is the ratio of number of white pixels on the edge to 
%   the number on white pixels in MASK. PENFLU is the fraction of 
%   fluorescence on the edge in masked positions.

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

pixelt = sum(protbin(:));
protclean = ml_mask(protimg, protbin);
fluort = sum(protclean(:));

pixele = 0;
fluore = 0;

for m = 1 : size(protimg, 3)
    bw = edge(protbin(:, :, m));
    pixele = pixele + sum(bw(:));
    tmp = ml_mask(protimg(:, :, m), bw);
    fluore = fluore + sum(tmp(:));
end

pencpixel = pixele / pixelt * 100;
penflu = fluore / fluort * 100;
