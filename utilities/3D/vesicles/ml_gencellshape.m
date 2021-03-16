function cellcode2 = ml_gencellshape(model,cellcode)
%ML_GENCELLSHAPE Generate a cell shape.
%   CELLCODE2 = ML_GENCELLSHAPE(MODEL,CELLCODE) returns a cell code (see 
%   ML_PARSECELL for the data structure.) which contains a generated cell
%   shape based on the the cell code CELLCODE. MODEL is a structure
%   returned from ML_RDISTPCA.

%   Ting Zhao
%
%   Copyright (c) 2006-2013 Murphy Lab
%   Carnegie Mellon University
%
%   May 13, 2013 I. Cao-Berg Updated method so that is has a maximum number
%   of resamples. This is useful for debugging and providing feedback to
%   the user
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

%maximum number of resamples if there is invalid and ratio and no minimal
%ratio
maximum_number_of_resamples = 10E6;

if nargin < 2
    error('Exactly 2 arguments are required');
end

isrepeat = 'y';
while isrepeat=='y'
    %Generate distance ratios.    
    x = ml_rnd(model.stat);
    
    %Check minimal ratios
    if ~isfield(model,'minratio')
         model.minratio = [];
    end
    
    %If the minimal ratio is used
    if ~isempty(model.minratio)
        if any(x<=model.minratio) %Set to minimal ratio
            if length(model.minratio==1)
                x(x<=model.minratio) = model.minratio;
            else
                x(x<=model.minratio) = model.minratio(x<=model.minratio);
            end
        end
    end
    
    %Try resampling if there is invalid ratio and no minimal ratio
    counter = 0;
    if isempty(model.minratio)
        if any(x<=1)
            counter = counter + 1;
            
            %icaoberg 13/5/2013
            if counter == maximum_number_of_resamples
                disp( 'Maximum number of resamples reached. Invalid ratios and no minimal ratio found.' );
                cellocde2 = {};
                return
            end
            
            continue;
        end
    end
    
    angles = 1:1:360/model.anglestep;

    %disp( 'Calculating the distances of cell boundary to nuclear center' );
    ddists = (x-1).*ml_shiftdist(cellcode.nucdist(angles), ...
        cellcode.nucmangle);
    ddists = circshift(ddists,[0 round(cellcode.nucmangle)]);
    
    %disp( 'Moving the center of nucleus to [0 0]' );
    nuchitpts = ml_addrow(...
        cellcode.nuchitpts(angles,:),-cellcode.nuccenter);
    
    %Deformation from a nucleus to a cell.
    pts2 = ml_moldshape(nuchitpts,ddists');
    
    cellcode2 = cellcode;
    cellcode2.da = model.anglestep;
    cellcode2.nuchitpts = nuchitpts;
    cellcode2.nucellhitpts = pts2;
    
    isrepeat = 'n'; %this feature is turned off because sometimes it is
                    %not supported by X-Win
end

