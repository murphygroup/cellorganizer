function seeds = fgseeds( I, PIXSIZE, DISKSIZE, HLEVELS, NUCSIZE,ARRATIO)
% SEEDS = FGSEEDS( I, PIXSIZE, DISKSIZE, HLEVELS, NUCSIZE)
%
% Inputs: I = the nucleus image
%         PIXSIZE = the number of microns per pixels (default 0.125)
%         DISKSIZE = the diameter in microns if the blurring kernel
%            (default 2)
%         HLEVELS = the local threshold tolerance (default 27)
%         NUCSIZE = range of acceptable nuclei sizes (default [10 40])
%
% Output: SEEDS = the seed image, to be used with getvoronoi.m or
%            seededwatershed.m
%
% Justin Newberg
%   for Murphy Lab
%

% Copyright (C) 2008  Murphy Lab
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

if ~exist( 'PIXSIZE','var')
    PIXSIZE = 0.125;
end
if ~exist( 'HLEVELS','var')
    HLEVELS = 27;
end
if ~exist( 'DISKSIZE','var')
    DISKSIZE = 2;
end
if ~exist( 'NUCSIZE','var')
    NUCSIZE = [10 40];
end
if ~exist('ARRATIO','var')
    ARRATIO = 0.95;
end

minnucsize = NUCSIZE(1);
maxnucsize = NUCSIZE(2);


f = fspecial('disk',DISKSIZE/PIXSIZE/2);
g = imfilter( I, f, 'replicate');

% seeds = ~imextendedmin(g, HLEVELS);
seeds = imextendedmax(g, HLEVELS);
seeds = imfill( seeds, 'holes');

bwl = bwlabel( seeds, 4);
stats = regionprops( bwl,'MajorAxisLength','Area');%,'ConvexArea');
majlen = zeros(length(stats),1);
% arratio = zeros(length(stats),1);
for i=1:length(stats)
    majlen(i) = stats(i).MajorAxisLength;
%    arratio(i)=stats(i).Area/stats(i).ConvexArea;
end
majlen = majlen*PIXSIZE;

seeds = 255*uint8(seeds);
%idx = ismember(bwl,find(arratio<ARRATIO));
%seeds(idx) = 32;
idx = ismember(bwl,find(majlen<minnucsize));
seeds(idx) = 64;
idx = ismember(bwl,find(majlen>maxnucsize));
seeds(idx) = 128;
idx = ismember(bwl,find(majlen<minnucsize/2));
seeds(idx) = 0;

idx = seeds & ~imclearborder(seeds>0);
seeds(idx) = 128;

return
