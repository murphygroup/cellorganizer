function filelisting = ml_dir( pattern )
% FILELISTING = ML_DIR( PATTERN) Returns file listing of wildcard PATTERN (e.g. '*.jpg').
% FILELISTING is a cell array of strings of filenames without path.
% the listing is sorted by name alphabetically

% Copyright (C) 2006-2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% 1-March-2013 Arun Sampath fixed a bug in sorting, use sort_nat function 
%                instead of natural sorting.
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

files = dir( pattern);
L = length( files );
if( L == 0)
  filelisting = {};
else
  [filelisting{1:L,1}] = deal( files.name);
  filelisting = sort_nat( filelisting );
end
