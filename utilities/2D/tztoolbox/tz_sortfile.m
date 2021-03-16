function files2 = tz_sortfile(files,pos)
%TZ_SORTFILE Sort the file names with their numbers.
%   FILES2 = TZ_SORTFILE(FILES) sorts the string array FILES in ascending 
%   order. The orders are decided by the number in each string in FILES.
%
%   FILES2 = TZ_SORTFILE(FILES,POS) specifies the position wheret the 
%   number is obtained from each string.
%
%   See ML_GETFILENUM for more details.
%   
%   See also ML_GETFILENUM

%   20-Jan-2006 Initial write T. Zhao
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

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('pos','var')
    pos = 1;
end

if isempty(files)
    files2 = {};
    return;
end

for i=1:length(files)
    fileNumbers(i) = ml_getfilenum(files{i},pos);
end

[tmp,sortIndices] = sort(fileNumbers);

files2 = files(sortIndices);
