function [ centrosome_coordinate ] = img2centrosome_coord( image )
%The centrosome is the brightest point in the MT image

%we are currently NOT implementing the method in Li, 2012:
%Centrosome location detection. The 3D coordinate of the centrosome was 
% estimated by breaking the problem into two parts. First, the XY-coordinate 
% was estimated and then the Z-coordinate. The XY-coordinate was chosen as 
% the pixel with the maximum intensity value in the vicinity of the nucleus 
% after smoothing with an averaging filter of size 25 pixels on the tubulin 
% channel image (as for cell image). For the Z-coordinate, we used linear 
% regression to estimate the location as a function of the following 
% predictor variables: (i) Maximum intensity of the microtubule image, (ii) 
% Mean intensity of the microtubule image, and (iii) pixel intensity of the 
% XY coordinate in the microtubule image. The coefficients of the linear 
% regression were estimated from the 3D HeLa images where the 3D centrosome 
% as described previously [8]. The estimated centrosome is then used to act 
% as an organizer for microtubules and all generated microtubules start from 
% it.

[~, ind] = max(image(:));
centrosome_coordinate = ind2sub_vec(size(image), ind);


end

