function pts = tz_img2crd(img,param)
%TZ_IMG2CRD Convert an image into coordinate shape.
%   PTS = TZ_IMG2CRD(IMG) returns a set of points that form a shape
%   extracted from the image IMG.
%   
%   PTS = TZ_IMG2CRD(IMG,PARAM) extracts a shape from IMG by a method
%   specified by PARAM. PARAM has following fields:
%       'method' - extraction method
%           'medaxis' - medial axis extraxtion
%           'contour' - contour tracing     
%   
%   See also

%   31-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = tz_initparam(param,struct('method','medaxis'));

switch param.method
    case 'medaxis'
        [imgaxis,axln,dists,borders] = tz_imaxis(img);
        pts1(:,1) = axln(:,1);
        pts1(:,2) = borders(:,1);
        pts2(:,1) = pts1(:,1);
        pts2(:,2) = borders(:,2);
        pts = [pts1;flipud(pts2)];
        pts(end+1,:) = pts(1,:);
    case 'contour'
        pts = tz_tracecontour(double(bwperim(img>0)));
end