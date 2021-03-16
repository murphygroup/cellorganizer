function img = tz_alignimgs(imgs,param)
%TZ_ALIGNIMGS Align multiple images to one image.
%   IMG = TZ_ALIGNIMGS(IMGS) returns an image that contains multiple images,
%   which are included in IMGS. IMGS is a cell array. Each of its elements
%   could be an [image] or the file path of an image. Currently it requires
%   all images have the same size and only have one channels.
%   
%   IMG = TZ_ALIGNIMGS(IMGS,PARAM) customizes the alignment by the structure
%   PARAM, which has the following fields:
%       'ncol' - number of columns (number of images in each row)
%       'width' - width of each image
%       'colspc' - column space
%       'rowspc' - row space
%       'res' - resolution
%       'label' - label of each image. default: [], which means no label.
%       'bg' - background value. Default: 0.
%
%   See also

%   14-Mar-2007 Initial write T. Zhao
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


if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param, ...
    struct('ncol',3,'width',100,'colspc',10,'rowspc',10,'res',[], ...
    'label',[],'bg',0,'paren','yes','ml_imaddtext',struct([]),'edge','no'));

img = [];
offset = [0 0];
nimg = length(imgs);
nrow = ceil(nimg/param.ncol);
% newimgsize = [param.width*param.ncol+param.colspc*(param.ncol-1), ...
%     param.*param.width+param.rowspc*(param.ncol-1)];

for i=1:length(imgs)
    if ischar(imgs{i})
        subimg = imread(imgs{i});
    else
        subimg = imgs{i};
    end
    
    imgsize = size(subimg);
    if ~isempty(param.width)
        scale = param.width/imgsize(2);
    else
        scale = 1;
    end
    
    if  strcmp(param.edge,'yes')
        labelimg = bwlabel(subimg>0);
        subimgs = {};
        nlabel = max(labelimg(:));
        if nlabel==1
            labelimg = bwlabel(bwmorph(subimg>0,'erode'));
            nlabel = max(labelimg(:));
        end
        
        for j=1:nlabel
            subimg2 = labelimg==j;
            subimgs{j} = ml_edgeresize(subimg2,scale);
            if j>=2 | nlabel<3
                subimgs{j} = imdilate(subimgs{j},ones(2,2));
            end
        end
%         subimgs{3} = imresize(labelimg==3,scale);
        subimg = zeros(size(subimgs{1}));
        for j=1:length(subimgs)
            subimg = subimg | subimgs{j};
        end
    else
        subimg = imresize(subimg,scale,'bicubic');
    end
    
    subimg = mat2gray(subimg);
    
    if ~isempty(param.label)
        if strcmp(param.paren,'yes')
            label = ['(' param.label ')'];
        else
            label = param.label;
        end
        subimg = ml_imaddtext(subimg,label,param.ml_imaddtext);
        param.label = char(param.label+1);
    end
    
    if i==1
        img = subimg;
    else
        img = tz_immerge(img,subimg,'min',offset,param.bg);
    end
    
    if mod(i,param.ncol)==0
        offset(1) = offset(1)+size(subimg,1)+param.rowspc;
        offset(2) = 0;
    else
        offset(2) = offset(2)+size(subimg,2)+param.colspc;
    end    
end

