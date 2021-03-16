function img2 = tz_imcropbg2(img,mask)
%TZ_IMCROPBG2 Crop part of background in a mask out
%   IMG2 = TZ_IMCROPBG2(IMG,MASK) returns an image by cropping background
%   out. Here background means pixels with value no greater than 0 in MASK.
%   
%   See also TZ_IMCROPBG

%   10-Nov-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

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

if nargin < 2
    error('Exactly 2 arguments are required')
end

img2=img;

img2(:,sum(mask,1)==0)=[];
img2(sum(mask,2)==0,:)=[];