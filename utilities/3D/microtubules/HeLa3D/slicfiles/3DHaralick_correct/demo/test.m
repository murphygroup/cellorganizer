addpath('../matlab');
addpath('../matlab/mex');

image(:,:,1) = imread('Synthetic_test1.bmp');
image(:,:,2) = imread('Synthetic_test2.bmp');
image (:,:,3) = imread('Synthetic_test3.bmp');
which ml_3Dcoocmat

direct = int32([1 1 0; 1 0 1; 0 1 1; 1 0 -1; 0 0 1]);
[f,n] = ml_3Dcoocmat(image,direct)
%[f,n] = ml_removegraylevel(f,n,0)
nb = int32(length(n))
[features,featuresnames] = ml_haralicktexture(single(f),n,nb);
features'
featuresnames
 
 