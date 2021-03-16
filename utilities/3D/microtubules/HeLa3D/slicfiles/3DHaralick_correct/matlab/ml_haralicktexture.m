function [features, featurenames] = ml_haralicktexture(coomat,grayLevels);
% COOCMAT = LM_3DCOOCMAT(IMAGE,direct,2,1);
%
% input: 
% - image: unsigned int8 image (gray-levels)
% - direction: matrix n-by-3 containing 0, 1 or -1. It gives the direction in x,y,z (respectively the first, second, third column).
%   n directions will be associated (merge) to build the cooccurrence matrix and the function will return only one cooccurrence matrix of size (height and width) equal to the gray levels in the image (zeros included).
%   This matrix has to be cast in int32 (eg. direction = int32(direction))
