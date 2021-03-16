function img2=tz_extendimg(img,k,g)

%function img2=tz_extendimg(img,k)

if ~exist('g','var')
    g=0;
end

imgsize=size(img);
img=[zeros(imgsize(1),k)+g,img,zeros(imgsize(1),k)+g];
imgsize=size(img);
img2=[zeros(k,imgsize(2))+g;img;zeros(k,imgsize(2))+g];