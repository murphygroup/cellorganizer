function ml_savefig(h,filepath,param)
%ML_SAVEFIG Save a figure into an image.
%   ML_SAVEFIG(H,FILEPATH) save the figure with handle H into the image
%   file FILEPATH.
%   
%   ML_SAVEFIG(H,FILEPATH,PARAM) specifies how to save the figure using the
%   structure PARAM, which has the following fields:
%       'size' - size of the image after saving
%       'paperSize' - printing paper size, only available for eps file
%       'scale' - scale the image intensity or not. not availabe for eps file.
%
%   See also

%   21-Jun-2006 Initial write T. Zhao
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
    error('2 or 3 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('size',[],'paperSize',[1275 1650], ...
                                  'mode','img','scale','yes'));

[pathstr,name,ext] = fileparts(filepath);

if strcmp(ext,'.fig')
    saveas(h,filepath,'fig');
    return;
end

if strcmp(ext,'.eps') | ~strcmp(param.mode,'img')
    if strcmp(ext,'.eps')
        param.mode = 'deps';
    end
    if ~isempty(param.size)
        normSize = fliplr(param.size)./param.paperSize;
        set(h,'PaperUnits','normalized');
        set(h,'PaperPosition',[0 0 normSize]);
    end
    print(h,['-' param.mode],filepath);
else
    if ~isempty(param.size)
        set(h,'Position',[20 20 fliplr(param.size)]);
    end
    X = getframe(h);
    
    if isempty(X.colormap)
        if strcmp(param.scale,'yes')
            X.cdata = mat2gray(X.cdata);
        end
        imwrite(X.cdata,filepath);
    else
        imwrite(X.cdata, X.colormap,filepath);
    end
end
