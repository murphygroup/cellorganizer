function [segdna, segcell] = ml_rescaleImage2Cell( segdna, segcell, stacknumcell )
%ML_RESCALEIMAGE2CELL

% Author: Ivan E. Cao-Berg (icaoberg@scs.cmu.edu)
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
ratio = stacknumcell/length(find(sum(sum(segcell))>=0));

%find new number of slices for both cell and dna images
%assumes segcell has already been trimmed to remove blank slices
newsize = round(ratio*size(segcell,3));

segdna = tp_stretch3d(segdna, newsize);
segcell = tp_stretch3d(segcell, newsize);
