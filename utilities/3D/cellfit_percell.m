function cellfit = cellfit_percell( seg_dna, seg_cell, param )
%Trains per-cell cell shape model

% Greg Johnson
%
% Copyright (C) 2016 Murphy Lab
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

%First check if they exist already

[cellcodes, equatorZ, cellnucheightratio, nucbottomslice] = find_cell_codes( seg_dna, seg_cell, param );

if ~isempty( cellcodes )
    [rad_ratio,nucdist,nuccelldist] = extract_radius_ratio(cellcodes, equatorZ);
    
    cellfit.cellcodes = cellcodes;
    cellfit.equatorZ = equatorZ;
    cellfit.cellnucheightratio = cellnucheightratio;
    cellfit.nucbottomslice = nucbottomslice;
    
    cellfit.rad_ratio = rad_ratio;
    cellfit.nucdist = nucdist;
    cellfit.nuccelldist = nuccelldist;
else
    if param.verbose
        disp( 'Cell codes is empty.' )
    end
    cellfit = [];
end
end