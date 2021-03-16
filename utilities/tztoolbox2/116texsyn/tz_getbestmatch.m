function [mwnd,mmask,mpos,mdist]=tz_getbestmatch(img,mask,wnd,qmask,df,maxmasked)

%function [wnd,pos,dist]=tz_getbestmatch(img,mask,wnd,df,minmasked)
%
%OVERVIEW:
%   Find the best matched texture window in an image
%PARAMETERS:
%   img - image for searching
%   mask - mask of the image
%   wnd - input query
%   df - distance function
%       'pixel'
%       'hist'
%   maxmasked - windows with number of masked pixels higher than maxmasked will be ignored
%RETURN:
%   mwnd - The best matched window
%   mmask - The mask on mwnd
%   mpos - position of the center of the window in the image
%   mdist - the minimal distance
%DESCRIPTION:
%

mdist = Inf;
mwnd = [];
mpos = [];

if isempty(mask)
    mask=ones(size(img));
end

k=floor(size(wnd,1)/2);

extimg=tz_extendimg(img,k,0);
extmask=tz_extendimg(mask,k,0);

for i=1:size(img,1)
    for j=1:size(img,2)
        
        maskwnd=tz_getnbwnd(extmask,[i,j],k,0);
        if ~isempty(qmask)
            maskwnd=maskwnd.*qmask;
        end
        
        maskwnd(k+1,k+1)=0;
        
        if sum(sum(maskwnd==0))<=maxmasked+1
            imgwnd=tz_getnbwnd(extimg,[i,j],k,0);
            tempdist=tz_wnddist(wnd,imgwnd,maskwnd,df);
            if tempdist<mdist
                mdist=tempdist;
                mwnd=imgwnd;
                mmask=maskwnd;
                mpos=[i,j];
            end
        end
    end
end

