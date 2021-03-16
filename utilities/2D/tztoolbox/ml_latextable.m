function ml_latextable(filename,data,names)
% ML_LATEXTABLE - print a matrix to a file for easy inclusion in latex
%

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

% M. Boland - 29 Apr 1999

% $Id: ml_latextable.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

if(~isstr(filename))
	error('FILENAME must be a string containing the full name of the file to be written') ;
end

if(ndims(data)>2)
	error('DATA must be a vector or 2D matrix') ;
end

if(length(names)>size(data,1) | ~iscell(names))
	error('NAMES must be a cell array and have the same number of elements as DATA has rows') ;
end

fid=fopen(filename,'w') ;
if(fid==-1)
	error('Unable to open the input file');
end

for i=1:size(data,1)
  fprintf(fid,'%s & ', names{i}) ;
  for j=1:size(data,2)
    fprintf(fid,'%0.0f%%',data(i,j)) ;
    if(j<size(data,2))
      fprintf(fid,' & ') ;
    else
      fprintf(fid,' \\\\[0.1in]\n') ;
    end
  end
end

fclose(fid) ;
