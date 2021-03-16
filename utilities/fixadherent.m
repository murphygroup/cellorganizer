function [rbot_slice, rtop_slice ] = fixadherent( segcell, bot_slice, top_slice )
%FIXADHERENT forces an adherent cell to have its largest slice on the bottom of the segmented cell
%
%Inputs:
% segcell = current segmented cell image
% bot_slice = bottom slice where signal exists
% top_slice = top slice where signal exists
%
%Outputs:
% rbot_slice = corrected bottom slice
% rtop_slice = corrected top slice

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

%Modified by:
% March 21, 2013 D. Sullivan, added license. 
% June 6, 2013 D. Sullivan, changed rbot_slice to match our current usage.

top = size(segcell,3);
areas = sum(sum(segcell));
[junk,truebottom] = max(areas);

%D. Sullivan 6/6/13, fixed rbot_slice problem
% I think this was correct when we pass only slices with fluorescence, but
% currently we pass all slices so the truebottom will have the correct
% indices
% rbot_slice = bot_slice + truebottom - 1;
rbot_slice = truebottom;

rtop_slice = top_slice;
end
