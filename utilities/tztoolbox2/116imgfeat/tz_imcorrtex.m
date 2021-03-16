function tex=tz_imcorrtex(img,mask,bin,option)
%TZ_IMCORRTEX Co-occurence matrix of an image. Obsolete.
%
%See also TZ_IMCORRTEX2

%   TEX = TZ_IMCORRTEX(IMG,MASK,BIN,OPTION) returns the co-occurence matrix
%   of the image IMG, which will be masked if MASK is not empty. BIN
%   specifies how many bins will be used for the matrix. OPTION is a cell
%   array of strings. If OPTION{1} is 'rm0', dark pixels (gray level 0)
%   will not be counted. If OPTION{2} is 'bin', IMG will be BINNED.

%   17-Sep-2005 Initial write T. Zhao
%   ??-???-???? Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%       - add comments
%       - change function name tz_texture --> tz_imcorrtex
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_imcorrtex','tz_imcorrtex2'));

if nargin < 4
    error('Exactly 4 arguments are required')
end

if ~isempty(mask)
   img=img.*mask+mask;
else
    if strcmp(option{1},'rm0')
        mask=img>0;
    else
        mask=ones(size(img));
    end
    
end

minlevel=min(img(mask==1));
maxlevel=max(img(mask==1));

if strcmp(option{2},'bin')
    binimg=round((img-minlevel)/(maxlevel-minlevel)*(bin-1)+1);
else
    binimg=img;
end

binimg=binimg.*mask;

tex=zeros(bin,bin);

imgsize=size(img);

for i=1:imgsize(1)
    index1=binimg(i,1:end-1);
    index2=binimg(i,2:end);
    for j=1:length(index1)
        if index1(j)~=0 & index2(j)~=0
            index=sort([index1(j),index2(j)]);
            tex(index(1),index(2))=tex(index(1),index(2))+1;
        end
    end
end

for i=1:imgsize(2)
    index1=binimg(1:end-1,i);
    index2=binimg(2:end,i);
    for j=1:length(index1)
        if index1(j)~=0 & index2(j)~=0
            index=sort([index1(j),index2(j)]);
            tex(index(1),index(2))=tex(index(1),index(2))+1;
        end
    end
end