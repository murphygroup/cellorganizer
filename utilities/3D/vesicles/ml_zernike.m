function [znames, zvalues] = ml_zernike(I,D,R)
% [ZNAMES, ZVALUES] = ML_ZERNIKE(I,D,R) Zernike moments through degree D 
% ML_ZERNIKE(I,D,R),
%     Returns a vector of Zernike moments through degree D for the
%     image I, and the names of those moments in cell array znames. 
%     R is used as the maximum radius for the Zernike polynomials.
%
%     For use as features, it is desirable to take the 
%     magnitude of the Zernike moments (i.e. abs(zvalues))
%
%     Reference: Teague, MR. (1980). Image Analysis via the General
%       Theory of Moments.  J. Opt. Soc. Am. 70(8):920-930.
%

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

% 19 Dec 98 - M.V. Boland
%

% $Id: ml_zernike.m,v 1.3 2006/06/27 13:33:47 tingz Exp $

znames = {} ;
zvalues = [] ;

%
% Find all non-zero pixel coordinates and values
%
[Y,X,P] = find(I) ;

%
% Normalize the coordinates to the center of mass and normalize
%  pixel distances using the maximum radius argument (R)
%
Xn = (X-ml_imgmoments(I,1,0)/ml_imgmoments(I,0,0))/R ;
Yn = (Y-ml_imgmoments(I,0,1)/ml_imgmoments(I,0,0))/R ;

%
% Find all pixels of distance <= 1.0 to center
%
k = find(sqrt(Xn.^2 + Yn.^2) <= 1.0) ;


for n=0:D,
    for l=0:n,
        if (mod(n-l,2)==0)
            znames = [znames cellstr(sprintf('Z_%i,%i', n, l))] ;
            zvalues = [zvalues ml_Znl(n, l, Xn(k), Yn(k), ...
                    double(P(k))/sum(P))] ;
        end
    end
end





