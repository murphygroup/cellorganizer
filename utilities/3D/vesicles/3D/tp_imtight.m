function [img2,bound] = tp_imtight(img)
% IMG2 = TP_IMTIGHT removes extra blank boundary of a 3D image

imgbw = img > 0;
nnzpixel = squeeze(sum(sum(imgbw)));
[maxfluo,supid] = max(nnzpixel);
supshape = imgbw(:,:,supid);
img = img(:,:,find(nnzpixel,1):find(nnzpixel,1,'last'));

rsum = sum(supshape,1);
L = find(rsum,1);
R = find(rsum,1,'last');
csum = sum(supshape,2);
T = find(csum,1);
B = find(csum,1,'last');

img2 = img(T:B,L:R,:);
bound = [T-1,B+1,L-1,R+1];