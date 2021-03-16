function answer = get_ometiff_dimensions( ometiff , param )
% function answer = get_ometiff_dimensions()

% Gets the dimensions of an OME.TIFF image, specified by the input param
%
% Ulani Qi (uhq@andrew.cmu.edu), Martha Cryan (mcryan@andrew.cmu.edu)
%
% Copyright (C) 2017 Murphy Lab
% Computational Biology Department
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.o
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

data = bfopen(ometiff); %opens image using bioformats
omeMeta = data{1,4}; %gets the metadata block
switch nargin
    case 1 % if no parameter specified, return all of them
        answer = [ omeMeta.getPixelsSizeX(0).getValue() ...
                   omeMeta.getPixelsSizeY(0).getValue() ...
                   omeMeta.getPixelsSizeZ(0).getValue() ...
                   omeMeta.getPixelsSizeC(0).getValue() ...
                   omeMeta.getPixelsSizeT(0).getValue() ];
    case 2 % parameter specified
        if strcmp(param, 'xyz') == 1
            answer = [ omeMeta.getPixelsSizeX(0).getValue() ...
                       omeMeta.getPixelsSizeY(0).getValue() ...
                       omeMeta.getPixelsSizeZ(0).getValue() ];
        elseif strcmp(param, 'c') == 1
            answer = [ omeMeta.getPixelsSizeC(0).getValue() ];
        elseif strcmp(param, 't') == 1
            answer = [ omeMeta.getPixelsSizeT(0).getValue() ];
        else
            warning('Invalid parameter: use xyz, c, t only');
        end
end

disp(answer);
end
