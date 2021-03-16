function savefile = ml_makelabelfig(img,param)
%ML_MAKELABELFIG Make a labeled figure.
%   ML_MAKELABELFIG(IMG,PARAM) put a label in the image
%   IMG and the the labeled image in a directory SAVEDIR. Where the file is
%   saved is decided by the structure PARAM, which has the following fields:
%       * 'paperid' - ID of the paper. It could be any string.
%       * 'figidx' - index of the figure. It could be any positive integer.
%       * 'label' - label of the figure, like 'a','b' etc.
%       * 'savedir' - the directory to save the figure
%       'imgformat' - format of the saved image. Default is 'tif'.
%       'paren' - add parenthesis or not. 'yes' (default) or 'no'.
%       'overwrite' - overwrite an existed file or not. 'yes' or 
%            'no' (default).
%       'crop' - crop area. 2x2 matrix. The 1st row is for topleft and the 2nd
%            row is for size. Default: empty (no crop). If the first two
%            row has 0 value, the topleft will be decided automatically.
%       'imgsize' - size of the labeled image. If it is not specified or empty,
%            the size will be the same as that of IMG. It could be 
%            [height width] or a number of resize ratio.
%       'ml_imaddtext' - how to add the label into the image. This is the
%            parameter structure for ML_IMADDTEXT. See ML_IMADDTEXT
%            for more details.
%   IMG could be an [image] or the file path of an image.
%
%   The saved file path is savedir/paperid.figidx.label.[imagformat], which
%   is the returned value of the function.
%
%   See also

%   27-Oct-2006 Initial write T. Zhao
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

param = ml_initparam(param,struct('imgsize',[],'overwrite','no', ...
        'paren','yes','imgformat','tif','ml_imaddtext',struct([]), ...
        'crop',[]));

separator = '.';
savefile = [param.savedir filesep param.paperid separator ...
            num2str(param.figidx) separator param.label '.' ...
            param.imgformat];

if exist(savefile,'file')
    if strcmp(param.overwrite,'yes')
        warning(['The file ' savefile ' is overwritten']);
    else
        warning(['THe file ' savefile ' already exists. Skip this ' ...
                            'file']);
        return;
    end
end

if ischar(img)
    img = imread(img);
end

if ~isempty(param.crop)
    img2 = zeros(param.crop(2,1),param.crop(2,2),size(img,3));
    img2 = eval([class(img) '(img2)']);
    if any(param.crop(1,:)<=0)
        topleft = [];
    else
        topleft = param.crop(1,:);
    end
    for i=1:size(img,3)
        [img2(:,:,i),rect] = ml_imcrop(img(:,:,i),struct('method','rect', ...
            'topleft',topleft,'size',param.crop(2,:)));
        if i==1
            topleft = rect(1,:);
        end
    end
    img = img2;
end

if ~isempty(param.imgsize)
    img = imresize(img,param.imgsize,'bilinear');
    img = mat2gray(img);
end

if strcmp(param.paren,'yes')
    img = ml_imaddtext(img,['(' param.label ')'],param.ml_imaddtext);
else
    img = ml_imaddtext(img,param.label,param.ml_imaddtext);
end

imwrite(img,savefile,param.imgformat);
