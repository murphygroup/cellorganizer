function [psf,resolution] = micro2psf(microscope)
% MICRO2PSF  Takes in a string name of a microscope for which to return a
% psf for, if 'none' is specified this will return an empty matrix 
%
% Input:
% microscope = string name for the microscope 
%
% Output:
% psf

% Author: Devin Sullivan (devins@andrew.cmu.edu)
% Edited: Ivan E. Cao-Berg (icaoberg@scs.cmu.edu)
% January 21, 2013 D. Sullivan added a resolution parameter to the psf file
% call
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

if nargin ~= 1
 error('micro2psf: Wrong number of input arguments');
end

if isempty( microscope )
 warning('micro2psf: No microscope specified.');
 psf = [];
end

try
switch lower(microscope)
    case 'svi'
        try
         load( 'data/HPA_285' );
        catch
         warning('micro2psf: Unable to load SVI psf');
         psf = [];
         resolution = [];
        end
    case 'none'
        psf = [];
        resolution = [];
    otherwise
        warning('micro2psf: Invalid microscope option.');
        psf = [];
        resolution = [];
end
catch
 warning('micro2psf: Invalid microscope option.');
 psf = [];
 resolution = [];
end
