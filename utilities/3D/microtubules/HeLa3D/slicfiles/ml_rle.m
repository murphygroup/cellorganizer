function enc_img = ml_rle( binimg)

% ENC_IMG = ML_RLE( BINIMG)
%
% Runlength encode a binary image of 0-s and 1-s.
% This results in better compression than, say PCX image format
% because binary input image is assumed, and therefore there is no
% need to record pixel values

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

tmp1 = uint8(binimg(1:end-1));
tmp2 = uint8(binimg(2:end));
change_pos = find( tmp1 ~= tmp2);
runlengths = uint32(diff( [0 change_pos]));
enc_img = struct('size',size(binimg),'runlengths',runlengths);
