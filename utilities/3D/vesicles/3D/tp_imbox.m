function bound = tp_imbox(img)
% IMG2 = TP_IMTIGHT removes extra blank boundary of a 3D image

% imgbw = img > 0;
% nnzpixel = squeeze(sum(sum(imgbw)));
% [maxfluo,supid] = max(nnzpixel);
% supshape = imgbw(:,:,supid);
% img = img(:,:,find(nnzpixel,1):find(nnzpixel,1,'last'));

%grj 4/23/15 removed the padding-by-one for images that abutt the
%boundaries of the image

warning('This function is depricated. Please use cropImgND.m.')

supshape = sum(img,3);
rsum = sum(supshape,1);
L = find(rsum,1);
R = find(rsum,1,'last');
csum = sum(supshape,2);
T = find(csum,1);
B = find(csum,1,'last');

%bound = [T-1,B+1,L-1,R+1];
bound = [T,B,L,R];
