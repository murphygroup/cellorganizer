function ml_makelibsvmfile(filepath,feats)
%ML_MAKELIBSVMFILE Make feature file for lib svm.
%   ML_MAKELIBSVMFILE(FILEPATH,FEATS) create a feature file with FILEPATH. 
%   This file contains the features from the [MCF] FEATS. The file has the
%   following format:
%       class feature_idx:value feature_idx:value ...
%       ...
%
%   See also

%   16-Jan-2007 Initial write S. Chen & T. Zhao
%   Copyright (c) 2007 Murphy Lab
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

nclass = length(feats);
nfeat = size(feats{1},2);

fid = fopen(filepath,'w');
for i=1:nclass
    for j=1:size(feats{i},1)
        fprintf(fid, '%d ',i);
        for k=1:nfeat
            fprintf(fid, '%d:%.6f ',k,feats{i}(j,k));
        end
        fprintf(fid, '\n');
    end
end
fclose(fid);
