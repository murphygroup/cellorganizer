function img2 = tz_nllapedge(img,wndsize)
%TZ_NLLAPEDGE Non-linear Laplace edge detection
%   IMG2 = TZ_NLLAPEDGE(IMG,WNDSIZE)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


%function img2 = tz_nllap_edge(img,wndsize)
%OVERVIEW:
%   non-linear Laplace edge detection
%PARAMETERS:
%   img - input image
%   wndsize - window size
%RETURN:
%   output image
%DESCRIPTION:
%   from the paper 'An Edge Detection Model Based on Non-linear Laplace Filtering'
%   Lucas J. et. al.
%
%HISTORY:
%   27-JUL-2003 Initial write TINGZ
%   14-DEC-2003 Modified TINGZ

[height width]=size(img);

if width==1
    loopj=1
else
    loopj=wndsize(2)+1:width-wndsize(2);
end

if height==1
    loopi=1;
else
    loopi=wndsize(1)+1:height-wndsize(1);
end

img2=img;
%img=tz_tophat(img,[0,3]);

for j=loopj
    for i=loopi
        maskwnd=img(i-wndsize(1):i+wndsize(1),j-wndsize(2):j+wndsize(2));
        img2(i,j)=max(maskwnd(:))+min(maskwnd(:))-2*img(i,j);
    end
end
