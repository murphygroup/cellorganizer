function binimg = ml_unrle( enc_img)

% BINIMG = ML_UNRLE( ENC_IMG)
%
% Unencode runlength encoded binary image of 0-s and 1-s.
% Image must have been coded with ML_RLE.

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

binimg = repmat(uint8(0),enc_img.size);
change_pos = cumsum(double(enc_img.runlengths));
startpos = change_pos(1:2:end)+1;
endpos = change_pos(2:2:end);
for i = 1:length(startpos)
    binimg(startpos(i):endpos(i)) = 1;
end
