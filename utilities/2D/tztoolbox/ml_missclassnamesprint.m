function mb_missclassnamesprint(filename,names,class)
% MB_MISSCLASSNAMESPRINT - prints NAMES next to the vector CLASS
%
%       Inputs:
%
%
%  M. Boland - 01 May 1999
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

% $Id: ml_missclassnamesprint.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

if (~iscell(names))
  error('NAMES must be a cell array') ;
end

if (length(names)<1)
  error('NAMES contains no elements') ;
end

if (length(names) ~= length(class))
  error('NAMES and CLASS must be of the same length') ;
end

fid=fopen(filename,'w') ;
if(fid==-1)
	error('Unable to open the output file');
end

for i=1:length(names)
	fprintf(fid,'%d\t%s\n',class(i),char(names(i))) ;
end

fclose(fid) ;

