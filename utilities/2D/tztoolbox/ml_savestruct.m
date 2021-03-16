function ml_savestruct(filename,s)
%ML_SAVESTRUCT Save structure in file
%   ML_SAVESTRUCT(FILENAME,S) save all fiels in S into the file with
%   file name FILENAME.

%   07-Apr-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

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
    error('Exactly 1 argumentis required')
end

c=struct2cell(s);
f=fieldnames(s);

savecmd=['save ' filename];

for i=1:length(c)
    eval([f{i} '=' 's.' f{i} ';']);
    savecmd=[savecmd ' ' f{i}];
end

eval(savecmd);

