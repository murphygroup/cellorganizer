function distr=tz_getdistgray(img)
%TZ_GETDISTGRAY Average gray levels of contours based on image edge.
%   DISTR = TZ_GETDISTGRAY(IMG) returns the average gray levels of pixels
%   on every contour which has the same distance to the edge. The edge is
%   detected from the binary version of IMG. DISTR has two columns with
%   each row for one contour. Each row has the form [DISTANCE,INTENSITY].
%   It is initially written for nucleus texture calibration

%   ??-???-???? Initial write T. Zhao
%   31-OCT-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

img2=tz_edgedist(img);

distvec=tz_mat2vec(img2);
imgvec=tz_mat2vec(img);

distvec(imgvec==0)=[];
imgvec(imgvec==0)=[];

distr=[];
while ~isempty(distvec)
    mindist=min(distvec);
    minpos=find(distvec==mindist);
    mingray=imgvec(minpos);
    distr=[distr;[mindist,mean(mingray)]];
    distvec(minpos)=[];
    imgvec(minpos)=[];
end
