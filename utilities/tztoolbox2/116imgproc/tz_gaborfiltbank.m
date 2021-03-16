function [out,avgout,ofreqs,othetas]=tz_gaborfiltbank(img,sigma)
%TZ_GABORFILTBANK Filter an image by gabor bank.
%   OUT = TZ_GABORFILTBANK(IMG,SIGMA)
%   
%   [OUT,AVGOUT,OFREQS,OTHETAS] = TZ_GABORFILTBANK(...)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


%function out=tz_gaborfiltbank(img,[])
%
%OVERVIEW:
%   filtering images by gabor bank
%PARAMETERS:
%   img - input image
%RETURN:
%   out - a series of filtered images
%   avgout - filtered images by combining images with the same freq
%   ofreqs - freq of each filter
%   othetas - orientation of each filter
%DESCRIPTION:
%   From paper Unsupervised Texture Segmentation Using Gabor Filters
%HISTORY:
%   09-NOV-2004 Initial write TINGZ
%   20-NOV-2004 Modified TINGZ
%       - add output avgout

bsize=max(size(img));
psize=2^ceil(log2(bsize));
% 
% if psize>size(img,2)
%     img=[img,zeros(size(img,1),psize-size(img,2))];
% end
% 
% if psize>size(img,1)
%     img=[img;zeros(psize-size(img,1),size(img,2))];
% end

thetas=[0,45,90,135];
freqs=2.^(0:log2(psize)-2)*sqrt(2)/psize;

nfilter=length(thetas)*length(freqs);
k=1;

sigmax=sigma;
sigmay=sigma;
maxout=0;
for i=1:length(freqs)
    avgout{i}=zeros(size(img));
    for j=1:length(thetas)
        [tmp,out{k}]=gaborfilter1(img,sigmax,sigmay,freqs(i),thetas(j));
        avgout{i}=avgout{i}+out{k};
        maxout=max(maxout,max(max(out{k})));
        ofreqs(k)=freqs(i);
        othetas(k)=thetas(j);
        k=k+1
    end
    avgout{i}=avgout{i}/length(thetas);
end

for i=1:length(out)
    out{i}=out{i}/maxout;
end

for i=1:length(avgout)
    avgout{i}=avgout{i}/maxout;
end