function img2=ml_pad(img)

%ML_PAD pad an image to reach the size 2^nx2^n, where n is the smallest power ...
% of two that is not less than the larger dimension of the image
%   ML_PAD(IMG)  appends zeros to the right and bottom of IMG

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

if s(1)<padsize
    img=[img;zeros(padsize-s(1),s(2))];
end

s=size(img);
if s(2)<padsize
    img=[img zeros(s(1),padsize-s(2))];
end

img2=img;