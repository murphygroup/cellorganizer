function rpixels=tz_genimgpixel(img,mask,n)

%function rpixels=tz_genimgpixel(img,mask,n)
%OVERVIEW:
%   Generate pixels from an image by its histogram
%PARAMETERS:
%   img - image for estimating pixel distributions
%   mask - mask for the image
%   n - number of pixels going to be generated
%RETURN:
%   rpixels - generated pixels

allpixels=reshape(img,1,prod(size(img)));

if ~isempty(mask)
    allmasks=reshape(img,1,prod(size(mask)));
    allpixels(allmasks==0)=[];
end

N=hist(allpixels,max(allpixels)-min(allpixels)+1);

p=N/sum(N);
[ignore,rpixels]=tz_mnornd(n,p,1);
rpixels=rpixels-1+min(allpixels);
