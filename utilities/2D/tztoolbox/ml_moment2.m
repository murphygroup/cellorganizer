function mom = ml_moment2(img)

%ML_MOMENT2 calculates central moments for an image
%   MOM=ML_MOMENT2(IMG) calculates moments for IMG
%   up to the sencond order. MOM is the structure.
%   The image center can be found in MOM.CX and MOM.CY.

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

%   24-Mar-2005 Initial write TINGZ



if min(img(:))==max(img(:))
    mom=struct('cx',size(img,1)/2,'cy',size(img,2)/2,'mu00',1,'mu11',0,'mu20', 0,'mu02',0);;
    warning('constant image');
    return;
end

[x,y]=find(img>0);
weights=img(find(img>0));
S=sum(weights);
centroid=sum([x.*weights,y.*weights],1)/S;
cpos=ml_addrow([x,y],-centroid);
% cpos=[x,y];

mu00=1;
mu11=sum(weights.*cpos(:,1).*cpos(:,2))/S;
mu20=sum(weights.*cpos(:,1).^2)/S;
mu02=sum(weights.*cpos(:,2).^2)/S;
% mu30=sum(weigths.*cpos(:,1).^3)/S;
% mu03=sum(weigths.*cpos(:,2).^3)/S;
mom = struct('cx',centroid(1),'cy',centroid(2),'mu00',mu00,'mu11',mu11,'mu20', mu20,'mu02',mu02);
