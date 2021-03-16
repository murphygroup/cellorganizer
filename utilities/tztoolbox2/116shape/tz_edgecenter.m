function [center,mangle] = tz_edgecenter(pts)
%TZ_EDGECENTER Obsolete. See ML_EDGECENTER.
%   CENTER = TZ_EDGECENTER(EDGE) returns the center of a solid object which 
%   has edge specified by a [curve] or an [image]  EDGE. If the number of
%   columns of EDGE is 2, EDGE will be take as a [curve], otherwise it will
%   be taken as an [image].
%   
%   [CENTER,MANGLE] = TZ_EDGECENTER(...) also returns the major angle of
%   the object. The unit is radian (same as the function TZ_BWMAJORANGLE).
%   
%   See also

%   24-Apr-2005 Initial write  T. Zhao
%   27-Jan-2006 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_edgecenter','ml_edgecenter'));

if nargin < 1
    error('Exactly 1 argument is required')
end

if size(pts,2)==2
    conpts=tz_showpts_2d(round(pts),'ln',1);
    edgeimg=tz_obj2img(conpts,[]);
else
    edgeimg = pts;
end

objimg=imfill(edgeimg,'hole');

% center=imfeature(objimg,'Centroid');
[x,y]=find(objimg==1);
center=mean([x,y],1);

mangle=tz_bwmajorangle(objimg);

% center(1)=sum(pts(:,1).^2)/sum(pts(:,1));
% center(2)=sum(pts(:,2).^2)/sum(pts(:,2));
