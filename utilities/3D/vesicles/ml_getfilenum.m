function num = ml_getfilenum(filename,pos)
%ML_GETFILENUM Extract number label from a file name.
%   ML_GETFILENUM(FILENAME) returns the number of the last occurrence from
%   FILENAME.
%
%   ML_GETFILENUM(FILENAME,POS) returns the number from FILENAME. The
%   occurrence of the number is specified by POS. If POS=1, it is the
%   last occurrence.
%
%   Example:
%       ml_getfilenum('image0403-1-023.tif') returns 23.
%       ml_getfilenum('image0403-1-023.tif',2) returns 1.
%       ml_getfilenum('image0403-1-023.tif',3) returns 403.

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

%   ??-???-???? Initial write T. Zhao
%   31-OCT-2004 Modified T. Zhao
%   15-MAR-2004 Modified T. Zhao
%       - add pos
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('At least one argument is required');
end

if nargin < 2
   pos = 1;
end

filename = fliplr(filename);
num = [];
curpos = 1;

for i = 1:length(filename)
   if filename(i) >= '0' & filename(i) <= '9'
      num = [num filename(i)];
   else
      if ~isempty(num)
         if curpos == pos
            break
         else
            curpos = curpos+1;
            num = [];
         end
      end
   end
end

if ~isempty(num)
   num = fliplr(num);
   num = str2num(num);
else
   warning('no number found');
   num = -1;
end
