function cellcode = ml_updatecellcode(cellcode,param)
%ML_UPDATECELLCODE Update cell code.
%   CELLCODE = ML_UPDATECELLCODE(CELLCODE,PARAM) returns a cell code that
%   is tranformed from the cell code CELLCODE based on the parameters
%   PARAM, which is a structure with the following fields:
%       'way' - how to update.
%           'cellhitpts' : cell hit points
%               'rdist' - distance ratio
%       
%
%   See also

%   03-Oct-2006 Initial write T. Zhao
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


if nargin < 2
    error('Exactly 2 arguments are required');
end

switch param.way
    case 'cellhitpts'
        ddists = (param.rdist-1).*ml_shiftdist(cellcode.nucdist, ...
            cellcode.nucmangle);
        ddists = circshift(ddists,[0 round(cellcode.nucmangle)]);

        %Move the center of nucleus to [0 0].
        cellcode.nuchitpts = ...
            ml_addrow(cellcode.nuchitpts,-cellcode.nuccenter);

        %Deformation from a nucleus to a cell.
        pts2 = ml_moldshape(cellcode.nuchitpts,ddists');

        cellcode.nucellhitpts = pts2;
        
        cellcode.nuchitpts = ...
            ml_addrow(cellcode.nuchitpts,cellcode.nuccenter);
        cellcode.nucellhitpts = ...
            ml_addrow(cellcode.nucellhitpts,cellcode.nuccenter);

    otherwise
        error(['Unrecognized way of updating: ' param.way]);
end
