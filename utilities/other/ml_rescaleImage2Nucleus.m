function [segdna, segcell] = ml_rescaleImage2Nucleus( segdna, segcell, stacknumdna )
%ML_RESCALEIMAGE2NUCLEUS

% Author: Robert F. Murphy (murphy@cmu.edu)
%
% Copyright (C) 2012 Murphy Lab
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


segdna = uint8(segdna);
segcell = uint8(segcell);

% find ratio of desired height to current number of non-blank slices in the dna image
ratio = stacknumdna/length(find(sum(sum(segdna))>=0));
%find new number of slices for both cell and dna images
%assumes segcell has already been trimmed to remove blank slices
newsize = round(ratio*size(segcell,3));

segdna = tp_stretch3d(segdna, newsize);
segcell = tp_stretch3d(segcell, newsize);

% get rid of any non-blank slices above and below the "real" nucleus due to round off
x=sum(sum(segdna));
for i=1:size(segdna,3)-stacknumdna+1
    sumx=sum(x(i:i+stacknumdna-1));
end

[junk,idx]=max(sumx);
idx=idx(1);
if idx>0 segdna(:,:,1:idx)=0; end
if idx+stacknumdna-1<=size(segdna,3) segdna(:,:,idx+stacknumdna-1:size(segdna,3))=0;
end
end
