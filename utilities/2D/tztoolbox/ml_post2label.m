function labels=ml_post2label(post)

%ML_POST2LABEL converts nx1 coding to contiguous integer labels
%   LABELS=ML_POST2LABEL(POST) converts each row of POST into an
%   integer, which is an element in POST. POST is a matrix for 
%   coding classes, e.g. class 1 is representd by [1 0 0 0 0] 
%   and class 3 is represented by [0 0 1 0 0] if there are five classes.
%   So POST is a nxk matrix if there are n samples and k classes 
%   and LABEL is a nx1 vector. 
%   
%   Example: If X = [1 0 0  ml_post2label(x) is [1
%                    1 0 0                       1
%                    0 1 0                       2
%                    0 0 1]                      3]
%
%   See also ML_LABEL2POST

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

%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   05-NOV-2004 Modified TINGZ
%       - add comments

if any(sum(post,2)~=1)
    error('wrong post')
end

[m,labels]=max(post,[],2);
% labels=labels';