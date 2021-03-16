function img2=ml_translate(img)

%ML_TRANSLATE translation of an image to its center
%   ML_TRANSLATE(IMG) pad IMG to a 2^nx2^n image. The mass center of the
%   image is at [2^(n-1),2^(n-1)]

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

s=size(img);
padsize=pow2(nextpow2(max(s)));
tcenter=padsize/2;

img2=img;

%calculate the mass center of the input image
totalmass=sum(img(:));
posarray=1:s(1);
xbar=round(posarray*sum(img,2)/totalmass);
posarray=1:s(2);
ybar=round(sum(img,1)*posarray'/totalmass);



%crop and extend
%top, left
if xbar>tcenter
    img2(1:xbar-tcenter,:)=[];
end

if ybar>tcenter
    img2(:,1:ybar-tcenter)=[];
end

if xbar<tcenter
    img2=[zeros(tcenter-xbar,size(img2,2));img2];
end

if ybar<tcenter
    img2=[zeros(size(img2,1),tcenter-ybar),img2];
end

%bottom, right
if size(img2,1)>padsize
    img2(padsize+1:end,:)=[];
end

if size(img2,2)>padsize
    img2(:,padsize+1:end)=[];
end

if size(img2,1)<padsize
    img2=[img2;zeros(padsize-size(img2,1),size(img2,2))];
end

if size(img2,2)<padsize
    img2=[img2,zeros(size(img2,1),padsize-size(img2,2))];
end



