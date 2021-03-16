function ml_writestringarray(varargin)
%ML_WRITESTRINGARRAY Write a string array into a file
%   ML_WRITESTRINGARRAY(OUTPUTFILE,S) writes the string array S into
%   OUTPUTFILE. Each element of S is a line. S could also be a vector.
%    
%   ML_WRITESTRINGARRAY(OUTPUTFILE,S1,S2,...) writes more than one string
%   array or vector into the file.

%   25-Jun-2005 Initial write TINGZ


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
    
if nargin < 2
    error('At least 2 arguments are required')
end

outputfile = varargin{1};

fid=fopen(outputfile,'w');
s = varargin{2};

for i=1:length(s)
    for j=2:nargin    
        if iscell(varargin{j})
            str = varargin{j}{i};
        else
            str = num2str(varargin{j}(i));
        end
        
        fprintf(fid,'%s',str);
        if j<nargin
            fprintf(fid,',');
        else
            fprintf(fid,'\n');
        end
    end
end

fclose(fid);
