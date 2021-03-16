function [img,ztight] = ic_imztight(img)
% IMG2 = TP_IMTIGHT removes extra blank boundary of a 3D image

imgbw = img > 0;
nnzpixel = squeeze(sum(sum(imgbw)));
[maxfluo,supid] = max(nnzpixel);
supshape = imgbw(:,:,supid);
img = img(:,:,find(nnzpixel,1):find(nnzpixel,1,'last'));
ztight = length(nnzpixel);