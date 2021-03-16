function [puncta_img,centers] = find_objects_local(img,flag,min_obj_size, ...
    max_obj_size, sigma, thresPerc )

% Copyright (C) 2020 Murphy Lab
% Computational Biology Department
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

disp(['Finding objects using local threshold']);
%     sigma = 5;
%     thresPerc = 0.1;
%     min_obj_size = 20; % for mitochondria
BW_pre=zeros(size(img));
puncta_img = BW_pre;
centers = {};
for i=1:size(img,3)
    %         background = imgaussfilt(img(:,:,i),sigma);
    BW_pre(:,:,i) = img(:,:,i)-imgaussfilt(img(:,:,i),sigma);
end
thres = thresPerc*max(max(BW_pre));
BW = BW_pre>thres;
objs = bwconncomp(BW);
if objs.NumObjects > 0
        centers_all = struct2cell(regionprops(objs,'Centroid'));
        centers_all = cellfun(@int16,centers_all,'UniformOutput',false);
        numPixels = cellfun(@numel,objs.PixelIdxList);
        idx=find(and(numPixels>min_obj_size,numPixels<max_obj_size));
        keep = zeros(size(idx));
        for j=1:length(idx)
            [xx,yy,zz]=ind2sub(size(img),objs.PixelIdxList{j});
            xwidth = max(xx)-min(xx);
            ywidth = max(yy)-min(yy);
            zwidth = max(zz)-min(zz);
            keep(j) = min([xwidth,ywidth,zwidth])>0;
        end
        idx = idx(find(keep));
        centers=centers_all(idx);
        for j=1:length(idx)
            puncta_img(objs.PixelIdxList{idx(j)}) = 1;
        end
end    
disp(['Object finding complete. Number of objects found is ' ...
        num2str(length(centers)) '.']);
end%find_objects_local
