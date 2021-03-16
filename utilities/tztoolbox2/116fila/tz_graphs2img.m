function img=tz_graphs2img(mts)
%TZ_GRAPHS2IMG Convert objects to an image.
%   IMG = TZ_GRAPHS2IMG(MTS) returns an image which contains all objects
%   in the one-level cell array MTS.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

allpts=[];
for i=1:length(mts)
    allpts=[allpts;mts{i}];
end

imgsize=round(max(allpts)-min(allpts));
imgsize=imgsize+[4,4];

center=round(imgsize/2);
allpts=round([allpts(:,1)+center(1),allpts(:,2)+center(2)]);

img=zeros(imgsize);
img(sub2ind(imgsize,allpts(:,1),allpts(:,2)))=1;