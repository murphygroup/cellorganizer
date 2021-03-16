function feats = ml_objtex3d(objects,param)
%ML_OBJTEX3D Texture features of a 3D object.
%   FEATS = ML_OBJTEX3D(OBJECTS) returns 26 texture features of the cell 
%   array of 3D objects OBJECTS.
%   
%   FEATS = ML_OBJTEX3D(OBJECTS,PARAM) specifies how to calculate the features 
%   by the structure PARAM, which has the following fields:
%       'orgres' - original resolution. Default: [] (no downsampling will be 
%            done).
%       'texres' - texture resolution. Default: [] (no downsampling will be 
%            done).
%       'bin' - number of bins for cooccurence matrix. Default: [] (original).
%   
%   See also

%   23-Mar-2007 Initial write T. Zhao
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

param = ml_initparam(param,struct('orgres',[],'texres',[],'bin',[]));
feats = [];

for i=1:length(objects)
    if isstruct(objects{i})
        obj = objects{i}.voxels;
    else
        obj = objects{i};
    end
    
    objimg = ml_obj2img(obj,[],{'3d','og'});
    
    if ~isempty(param.orgres) & ~isempty(param.texres)
        ratio = param.texres./param.orgres;
    else
        ratio = [];
    end
    
    err = 0;
    if ~isempty(ratio)
        if (length(find(ratio - round(ratio))))
            objimg = ml_3dimresize(objimg, 1/ratio(1), 1/ratio(3));
        else
            objimg = ml_downsize(objimg, ratio,text.method);
        end
    end

    if ~isempty(param.bin)
        objimg = uint8(floor(double (objimg) * (param.bin - 1) /... 
            double(max(objimg(:)))));
    end
    if size(objimg,1)==1 | size(objimg,2)==1
        err = 1;
    else
        eval('z = ml_3Dtexture(objimg);', 'err = 1;')
    end
    
    if (err)
        warning('problem with calculating features.');
        feat = zeros(1,26)+NaN;
    else 
        feat = [z(1:13, 14)' z(1:13, 15)'];
    end
    feats = [feats; feat];
end
