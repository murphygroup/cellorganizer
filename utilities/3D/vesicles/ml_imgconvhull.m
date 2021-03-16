function hull = ml_imgconvhull(bwimage)
% ML_IMGCONVHULL returns the convex hull of BWIMAGE as an image
% ML_IMGCONVHULL(IMAGE), where BWIMAGE is the binary image 
%    to be processed.
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

% 07 Aug 98 - M.V. Boland

% $Id: ml_imgconvhull.m,v 1.3 2006/06/27 13:33:47 tingz Exp $

if (isempty(bwimage))
    error('Invalid input image') ;
end

%
% Find all non-zero points in bwimage, and then identify the points on 
%   the convex hull
%
[Ty, Tx] = find(bwimage) ;
if (size(Tx,1)>2 & size(Ty,1)>2)
    succ=1;
	eval('k = convhull(Tx,Ty);','succ=0;');
    if succ==0
        k=1:length(Tx);
    end
    
    %Actually it's not suitable to use roipoly here
	hull = roipoly(bwimage, Tx(k), Ty(k)) ; 
else
    if length(Tx)==1
        hull=bwimage;
    else    
        pts=ml_getlinepts([Ty(1) Tx(1)],[Ty(2) Tx(2)]);
	    hull = ml_setimgptspixel(bwimage,pts);%zeros(size(bwimage,1), size(bwimage,2));
    end
end

%
% Convert the convex hull polygon to an image
%