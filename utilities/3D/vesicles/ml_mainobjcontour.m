function img2=ml_mainobjcontour(img)
%ML_MAINOBJCONTOUR Find the boundary of the biggest object in an image.
%   IMG2 = ML_MAINOBJCONTOUR(IMG) returns an image that contains the
%   boundary of the biggest object in the image IMG.
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% Copyright (C) 2007  Murphy Lab
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

if isempty(img)
    img2=[];
    return;
end

img=double(img>0);

img2=imfill(img,'holes');
imgedge=bwperim(img2,4);

img2=double(img2)-double(imgedge);

obj=ml_findmainobj2d_bw(img2,4);
%obj=ml_findmainobj_bw(img2,4);

img2=ml_obj2img(obj,size(img));
img2=img2==0;
img2=[ones(1,size(img2,2));img2;ones(1,size(img2,2))];
img2=[ones(size(img2,1),1),img2,ones(size(img2,1),1)];
img2=bwperim(img2,4);
img2(:,1)=[];
img2(:,end)=[];
img2(1,:)=[];
img2(end,:)=[];
