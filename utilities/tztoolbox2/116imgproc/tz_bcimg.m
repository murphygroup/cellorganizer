function img2=tz_bcimg(img,grayin,grayout)
%TZ_BCIMG Obsolete. See ML_BCIMG.
%   IMG2 = TZ_BCIMG(IMG,[MINV,MAXV],[LOW,HIGH]) rescales IMG to the range
%   [LOW HIGH], which is mapped from [MINV,MAXV].
%   IMG2 = TZ_BCIMG(IMG,[MINV,MAXV],[]) maps the minimal value and maximum
%   value in IMG to [LOW,HIGH].

%   27-SEP-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_bcimg','ml_bcimg'));

if nargin < 3
    error('Exactly 3 arguments are required')
end

img=double(img);

if isempty(grayin)
    grayin=[min(img(:)),max(img(:))];
end

if isempty(grayout)
    grayout=[0 255];
end

img2=(img-grayin(1))/(grayin(2)-grayin(1))*(grayout(2)-grayout(1))+grayout(1);
img2(img<=grayin(1))=grayout(1);
img2(img>=grayin(2))=grayout(2);
