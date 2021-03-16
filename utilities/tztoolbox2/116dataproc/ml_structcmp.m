function v = ml_structcmp(struct1,struct2)
%ML_STRUCTCMP Compare two structures.
%   V = ML_STRUCTCMP(STRUCT1,STRUCT2) returns 1 if the two structures
%   STRUCT1 and STRUCT2 are exactly the same. Otherwise it returns 0.
%   
%   See also

%   07-Nov-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 2
    error('Exactly 2 arguments are required');
end

fields1 = sort(fieldnames(struct1));
fields2 = sort(fieldnames(struct2));

if length(fields1)~=length(fields2)
    v = 0;
    return;
end

for i=1:length(fields1)
    if strcmp(fields1{i},fields2{i})~=1
        v = 0;
        return;
    end
    v1 = getfield(struct1,fields1{i});
    v2 = getfield(struct2,fields1{i});
    
    if istruct(v1) & istruct(v2)
        v = ml_structcmp(v1,v2);
        if v==0
            return
        end
    else
        if isstruct(v1) | istruct(v2)
            v = 0;
            return;
        end
    end
    
    if length(v1(:))~=length(v2(:))
        v = 0;
        return;
    end
    
    if any(v1(:)~=v2(:))
        v=0;
        return;
    end
    
end


v = 1;
