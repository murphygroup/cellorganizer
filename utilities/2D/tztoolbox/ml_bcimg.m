function img2=ml_bcimg(img,grayin,grayout)
%ML_BCIMG Change the brightness and constrast of an image.
%   IMG2 = ML_BCIMG(IMG,[MINV,MAXV],[LOW,HIGH]) rescales IMG to the range
%   [LOW HIGH], which is mapped from [MINV,MAXV].
%   IMG2 = ML_BCIMG(IMG,[MINV,MAXV],[]) maps the minimal value and maximum
%   value in IMG to [LOW,HIGH].

%   27-SEP-2004 Initial write T. Zhao
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

if nargin < 3
    error('Exactly 3 arguments are required')
end

img=double(img);

if isempty(grayin)
    grayin=[min(img(:)),max(img(:))];
end
grayin = double(grayin);
if isempty(grayout)
    grayout=[0 255];
end
grayout = double(grayout);

if sum(grayin~=grayout)>0
    img2=(img-grayin(1))/(grayin(2)-grayin(1))*(grayout(2)-grayout(1))+grayout(1);
    img2(img<=grayin(1))=grayout(1);
    img2(img>=grayin(2))=grayout(2);
else
    img2 = img;
end