function savedatfile(savepath,resolution,downsample)
%This function saves .dat file for loading primitives into blender
%currently it saves a single .dat file for each endosome with only size
%data.
%
%Inputs:
%savepath = string pointing to the desired folder to save the .dat files
%
%Outputs:
%set of .dat files

% Author: Devin Sullivan
%
% Copyright (C) 2013-2014 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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

mkdir(savepath)

if isempty(resolution)
    resolution = [1,1,1];
end

if isempty(downsample)
    downsample = [1,1,1];
end

if size(downsample)==1
    downsample = [downsample,downsample,downsample];
end

%first load temp results
load([pwd filesep 'temp' filesep 'primitives']);
% objsizevec = objsizevec.*repmat(resolution,[size(objsizevec,1),1]);
% objsizevec = objsizevec.*repmat(downsample,[size(objsizevec,1),1]);

filename = [savepath '.dat'];

for i = 1:size(objsizevec,1)
    dlmwrite(filename,['Endosome[',num2str(i),']'],'delimiter','','-append');
    dlmwrite(filename,objsizevec(i,:),'delimiter',' ','-append');
    dlmwrite(filename,objposvec(i,:),'delimiter',' ','-append');
    dlmwrite(filename,objrotvec(i,:),'delimiter',' ','-append');
    %     dlmwrite(filename,'cube','delimiter','','-append');
    %     dlmwrite(filename,objsizevec(i,:)*2,'delimiter',' ','-append');
end