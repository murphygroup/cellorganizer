function [mixes, objintens, offsets] = compile_gaussobjs(protparam)
% Fit gaussian distributions to the objects

% Tao Peng
%
% Copyright (C) 2012-2013 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% April 26, 2012 I. Cao-Berg Added new flag to save method to allow saving of
%              large matrices
% July 30, 2012 I. Cao-Berg Added a new parameter minObjSize for filtering 
%              objects by size
% August 1, 2012 R.F.Murphy Add debug code; remove unnecessary ct indix
% Jan 30, 2012 I. Cao-Berg Made declaration of temporary folder platform
%              and operating system independent
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
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



    %Check for temp results directory with per-cell results
   

%Now make a cell array of files 
files = ml_dir([temporaryFilesDirectory filesep 'gaussobjs_*.mat']);
mixes_tot = cell(1,length(files));
objintens_tot = cell(1,length(files));
offsets_tot = cell(1,length(files));

for i = 1:1:length(files)
    %load(['./temp/protein_objects_gaussian/original_objects/' protype '/obj' num2str(i) '.mat'])
    load([temporaryFilesDirectory filesep files{i}]);
    disp(['file: ' files{i}]) 

        %D. Sullivan 6/12/13 - now just compiling from the read temp
        %files
        mixes_tot{i} = mixes{1};
        objintens_tot{i}= objintens{1};
        offsets_tot{i}= offsets{1};
%             mixes{end+1} = cellmixes;
%             objintens{end+1}= cellobjintensity;
%             offsets{end+1}= cellobjoffset;

end 
mixes = mixes_tot;
objintens = objintens_tot;
offsets = offsets_tot;

