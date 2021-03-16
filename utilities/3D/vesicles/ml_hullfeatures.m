function [names, values, slfnames] = ml_hullfeatures(imageproc, imagehull)
% MB_HULLFEATURES - calculates features using the convex hull of IMAGE
% [NAMES, VALUES, SLFNAMES] = MB_HULLFEATURES(IMAGEPROC, IMAGEHULL) 
%     IMAGEPROC is the processed fluorescence image, IMAGEHULL is its binary 
%     convex hull.
%
% 12 Aug 98 - M.V. Boland
% June 2, 2002 - M.Velliste: added SLF feature names

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

% $Id: ml_hullfeatures.m,v 1.3 2006/06/27 13:33:47 tingz Exp $


% if ~isbw(imagehull)
%     error('The convex hull image must be binary.') ;
% end

%
% Initialize the return variables
%
names = {} ;
slfnames = {} ;
values = [] ;

%
% Fraction of the convex hull occupied by fluorescence
%    Ahull = area of the convex hull.
%
Ahull = bwarea(imagehull) ;
hullfract = length(find(imageproc))/Ahull ;
names = [names cellstr('convex_hull:fraction_of_overlap')] ;
slfnames = [slfnames cellstr('SLF1.14')] ;
values = [values hullfract] ;

%
% 'Shape factor' of the convex hull.
%    Phull = approximation to the hull perimeter.
%
Phull = bwarea(bwperim(imagehull)) ;
hullshape = (Phull^2)/(4*pi*Ahull) ;
names = [names cellstr('convex_hull:shape_factor')] ;
slfnames = [slfnames cellstr('SLF1.15')] ;
values = [values hullshape] ;

%
% Central moments of the convex hull
%
hull_mu00 = ml_imgcentmoments(imagehull,0,0) ;
hull_mu11 = ml_imgcentmoments(imagehull,1,1) ;
hull_mu02 = ml_imgcentmoments(imagehull,0,2) ;
hull_mu20 = ml_imgcentmoments(imagehull,2,0) ;

%
% Parameters of the 'image ellipse'
%   (the constant intensity ellipse with the same mass and
%   second order moments as the original image.)
%   From Prokop, RJ, and Reeves, AP.  1992. CVGIP: Graphical
%   Models and Image Processing 54(5):438-460
%

hull_semimajor = sqrt((2 * (hull_mu20 + hull_mu02 + ...
    sqrt((hull_mu20 - hull_mu02)^2 + ...
    4 * hull_mu11^2)))/hull_mu00) ;

hull_semiminor = sqrt((2 * (hull_mu20 + hull_mu02 - ...
    sqrt((hull_mu20 - hull_mu02)^2 + ...
    4 * hull_mu11^2)))/hull_mu00) ;

hull_eccentricity = sqrt(hull_semimajor^2 - hull_semiminor^2) / ...
    hull_semimajor ;

names = [names cellstr('convex_hull:eccentricity')] ;
slfnames = [slfnames cellstr('SLF1.16')] ;

values = [values hull_eccentricity] ;

