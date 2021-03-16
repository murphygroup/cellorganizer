function obj=tz_findmainobj_bw(img,n)
%TZ_FINDMAINOBJ_BW Obsolete. See ML_FINDMAINOBJ_BW.
%   OBJ = TZ_FINDMAINOBJ_BW(IMG) returns a 2-column matrix representing 
%   the extracted object from a binary image IMG. 
%   OBJ = TZ_FINDMAINOBJ_BW(IMG,N) specifies the parameter for BWLABEL. 
%   N defaults to 8.

%   30-SEP-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_findmainobj_bw','ml_findmainobj_bw'));

if nargin < 1
    error('1 or 2 arguments are required')
end

if nargin<2
    n=8;
end

limg=bwlabel(img,n);
objnum=max(limg(:));
        
lhist=[];
for i=1:objnum
    lhist(i)=sum(sum(limg==i));
end
[y,maxl]=max(lhist);
img(limg~=maxl)=0;

[r,c]=find(img==1);
obj=[r,c];