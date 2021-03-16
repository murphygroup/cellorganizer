function [E,mse,ssd,numpixels] = tz_errorsurface(T,rows,cols,J,I)
%TZ_ERRORSURFACE Modified from ERRORSURFACE by T. Zhao
%
% |----------------------------------------------------------|
% | Hybrid Texture Synthesis MATLAB package                  |
% |                                                          |
% | Author: Andrew Nealen                                    |
% |         Discrete Geometric Modeling Group                |
% |         Technische Universitaet Darmstadt, Germany       |
% |                                                          |
% | Note:   This is part of the prototype implementation     |
% |         accompanying our paper/my thesis                 |
% |                                                          |
% |         Hybrid Texture Synthesis. A. Nealen and M. Alexa |
% |         Eurographics Symposium on Rendering 2003         |
% |                                                          |
% |         Hybrid Texture Synthesis. A. Nealen              |
% |         Diplomarbeit (MSc), TU Darmstadt, 2003           |
% |                                                          |
% |         See the paper/thesis for further details.        |
% |----------------------------------------------------------|
%
% File errorsurface.m
%   This subroutine returns the error surface between a source (I) and
%   target (T) texture patch
%
%   [E,mse,ssd,numpixels] = errorsurface(T,rows,cols,J,I)
%
%   INPUT:
%   T  - the input texture
%   rows - patch size in row pixels
%   cols - patch size in column pixels
%   J - binary image mask support
%   I - image mask
%
%   OUTPUT:
%   E - the error surface (values in [0,1])
%   mse - the mean square difference (in [0,1])
%   ssd - the summed square difference
%   numpixels - number of pixels in the overlap region
%

% symbolic constants
RED = 1;
GREEN = 2;
BLUE = 3;
RED_WEIGHT = 0.299;
GREEN_WEIGHT = 0.587;
BLUE_WEIGHT = 0.114;

% compute the error surface E within the mask support J
E = zeros(rows,cols);
numpixels = 0;
ssd = 0;

numpixels = sum(J(:));
colorWeights = [RED_WEIGHT GREEN_WEIGHT BLUE_WEIGHT];
I2 = I(1:rows,1:cols,:);
J2 = J(1:rows,1:cols);
T2 = T(1:rows,1:cols,:);
E = RED_WEIGHT*J2.*(I2(:,:,RED)-T2(:,:,RED)).^2+ ...
    GREEN_WEIGHT*J2.*(I2(:,:,GREEN)-T2(:,:,GREEN)).^2+ ...
    BLUE_WEIGHT*J2.*(I2(:,:,BLUE)-T2(:,:,BLUE)).^2;

ssd = sum(E(:));

% build error surface E = (I - shiftT)^2, weighted by color channel weights.
% this supports the human visual system AND also maps our error to [0,1] when
% r,g,b in [0,1]
% for jj=1:rows,
%     for ii=1:cols,
%         % if we have mask support here, compute error surface E for this pixel
%         if (J(jj,ii) == 1),
%             numpixels = numpixels + 1;
% 
%             colorerror = RED_WEIGHT * (I(jj,ii,RED) - T(jj,ii,RED))^2;
%             E(jj,ii) = E(jj,ii) + colorerror;
%             ssd = ssd + colorerror;
%             
%             colorerror = GREEN_WEIGHT * (I(jj,ii,GREEN) - T(jj,ii,GREEN))^2;
%             E(jj,ii) = E(jj,ii) + colorerror;
%             ssd = ssd + colorerror;
%             
%             colorerror = BLUE_WEIGHT * (I(jj,ii,BLUE) - T(jj,ii,BLUE))^2;
%             E(jj,ii) = E(jj,ii) + colorerror;
%             ssd = ssd + colorerror;
%             
%         end
%     end
% end

if (numpixels < 1), mse = 0; else mse = ssd/numpixels; end
