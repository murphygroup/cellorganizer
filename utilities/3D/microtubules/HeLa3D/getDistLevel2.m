function [levels] = getDistLevel2(image,mask,numlevels)
% GETDISTLEVEL2

% Copyright (C) 2015-2016 Murphy Lab
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

% if ~exist('check','var')
%    check = 1;
% end

% [protim3,Dbothfin,segdna,segcell,dnaim3,cellim3,imgcent_coordinate] = getrealimage_hela(mask);
% if check && ismember(mask,[1,6,14,19,31,52])  %%
%    segcell = cat(3,segcell,zeros(size(segcell,1),size(segcell,2),16-size(segcell,3)));
% end
if isa( image, 'uint8' )
    mask = uint8(mask);
end

image = mask.*image;

[D,L] = bwdist(~mask);
unique_D = unique(D);
unique_D(1) = [];

if ~exist('numlevels','var')
    numlevels = 10;  %% Please also check hist_DistLevel_Curvelet_DistLevel.m for the compatibility.
    numlevels = min(numlevels,length(unique_D));
    %numlevels = max(numlevels,ceil(sqrt(length(unique_D))));
end
interval = floor(length(unique_D)/numlevels);

levels = cell(numlevels,1);

for i = 1:numlevels-1
    [idxr,idxc,idxz] = ind2sub(size(D),find(D>unique_D(1+interval*(i-1))-1e-6 & D<unique_D(interval*i+1)-1e-6));
    idx = find(D>unique_D(1+interval*(i-1))-1e-6 & D<unique_D(interval*i+1)-1e-6);
    ind = find(image(idx)); idxr = idxr(ind); idxc = idxc(ind); idxz = idxz(ind); idx = idx(ind);  %%The difference between getDistLevel2 and getDistLevel.
    levels{i} = [idxr(:),idxc(:),idxz(:),image(idx)];
    if isempty(levels{i})  %%
        warning('Empty level!');
        levels{i} = [0,0,0,0];
    end
end

i = numlevels;
[idxr,idxc,idxz] = ind2sub(size(D),find(D>unique_D(1+interval*(i-1))-1e-6 & D<unique_D(end)+1e-6));
idx = find(D>unique_D(1+interval*(i-1))-1e-6 & D<unique_D(end)+1e-6);
ind = find(image(idx)); idxr = idxr(ind); idxc = idxc(ind); idxz = idxz(ind); idx = idx(ind);  %%The difference between getDistLevel2 and getDistLevel.
levels{i} = [idxr(:),idxc(:),idxz(:),image(idx)];
if isempty(levels{i})  %%
    warning('Empty level!');
    levels{i} = [0,0,0,0];
end
