function seeds = bgseeds( I, RADIUS, HLEVELS, PIXSIZE)

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


if ~exist( 'HLEVELS','var')
    HLEVELS = 3;
end
if ~exist( 'RADIUS','var')
    RADIUS = 20;
end
if ~exist('PIXSIZE','var')
    PIXSIZE = 0.25;
end

minsize = 150*0.25;

f = fspecial('disk',RADIUS);
g = imfilter( I, f, 'replicate');

seeds = imextendedmin(g, HLEVELS);
bwl = bwlabel(seeds,4);

stats = regionprops(bwl,'MajorAxisLength');

majlen = zeros(length(stats),1);
for i=1:length(stats)
    majlen(i) = stats(i).MajorAxisLength;
end
majlen = majlen*PIXSIZE;


idx = ismember(bwl,find(majlen<minsize));
seeds(idx) = 0;

seeds = uint8(seeds);

return
