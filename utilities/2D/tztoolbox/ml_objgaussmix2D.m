function mix = ml_objgaussmix(obj,param)
%ML_OBJGAUSSMIX Fit gaussian mixture for an object.
%   MIX = ML_OBJGAUSSMIX(OBJ) returns a gaussian mixture structure that can
%   fit the [object] OBJ.
%   
%   MIX = ML_OBJGAUSSMIX(OBJ,PARAM) specifys how to fit the object by
%   PARAM:
%       'filter' - filtering the object. It is empty if there is no
%           filtering.
%       'mindist' - specify the minimal distance between the initialized
%           centers.
%       'isshow' - show the results (1) or not (0)
%       'gmm' - a structure to specify gmm parameters (see GMM):
%           'covartype' - covariance type
%           'ppca_dim' - pca dimension
%           'options' - options for specifying EM algorithm (see ML_GMMEM)
%
%   See also ML_OBJGAUSS

%   11-Jul-2006 Initial write T. Zhao
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

if nargin < 1
    error('1 or 2 arguments are required');
end

if size(obj,1)==1
    pts = obj(:,1:end-1);
    mix = gmm(size(pts,2),1,'full',1);
    mix.ncentres = 1;
    mix.centres = pts;
    mix.covars = zeros(size(pts,2));
    disp('The object has only one pixel');
    return;
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('filter',[],'mindist',0,'isshow',0, ...
    'gmm',struct([])));
param.gmm = ml_initparam(param.gmm, ...
    struct('covartype','full','ppca_dim',1,'options',zeros(1,14)));

[objimg,offset] = ml_obj2img2D(obj,[]);

if ~isempty(param.filter)
    dimobjimg = imfilter(objimg,param.filter,'same');
else
    dimobjimg = objimg;
end

% newidx = sub2ind(size(objimg),obj(:,1)+offset(1),obj(:,2)+offset(2));
% obj(:,3) = dimobjimg(newidx);

pts = ml_imlocalmax(dimobjimg);
if param.mindist>0
    pts = ml_mergepts(pts,param.mindist);
end

% if param.isshow==1
%     subplot(1,2,1)
%     imshow(objimg,[])
%     hold on
%     if ~isempty(pts)
%         plot(pts(:,2),pts(:,1),'x');
%     end
%     
%     hold off
%     drawnow
% end

scale = 100;

if ~isempty(pts)
    mix = gmm(size(pts,2),size(pts,1),param.gmm.covartype, ...
        param.gmm.ppca_dim);
    mix.centres = pts/scale;
else
    mix = gmm(size(obj,2)-1,1,param.gmm.covartype, ...
        param.gmm.ppca_dim);
end

[mix, options, errlog] = ml_gmmem2D(mix,obj(:,1:2)/100, ...
    param.gmm.options,obj(:,3));

mix.centres = mix.centres*scale;
mix.covars = mix.covars*scale^2;

if strcmp(mix.covar_type,'spherical')
    mix.covar_type = 'full';
    for i=1:mix.ncentres
        covars(:,:,i) = diag([mix.covars(i) mix.covars(i)]);
    end
    mix.covars = covars;
end

for i=1:size(mix.covars,3)
    zerovaridx = find(diag(mix.covars(:,:,i))==0);
    for j=1:length(zerovaridx)
        mix.covars(zerovaridx,zerovaridx,i) = 1e-2;
    end
end

% if param.isshow==1
%     img = ml_gaussmiximg(mix);
%     subplot(1,2,2)
%     imshow(img,[]);
% end
