function [out] = errorimage_src(data,I,J)
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
% File errorimage_src.m
%
%   [out] = errorimage_src(data,I,J)
%
%   This subroutine will compute the error
%   surface for a given texure (stored in 'data'), image mask I and mask 
%   support function J. 
%
%   extended to SOURCE IMPORTANCE WEIGHTING using a 
%   (source) texture importance map S from hybrid texture synthesis thesis
%
%   'Hybrid Texture Synthesis', Nealen A.
%   Diplomarbeit (MSc thesis), Technische Universitaet Darmstadt
%
%   complexity O(N log N) with N = w x h
%
%   INPUT:
%   data - the input datastructure as defined in hybridsynthesize.m and gen_input.m
%   I - the image mask
%   J - the binary mask support function
%
%   OUTPUT:
%   out - the error image with error values in [0,1]
%
%   The routine will return the error surface. indices
%   are 1-based, meaning the value in the upper left corner
%   is the error for shifting mask by +1 pixel in x and y 
%   directions. the zero value for no translation is stored
%   at index [h,w] (bottom right), which is eqivalent to no 
%   translation as the input texture is handled toroidally
%
%   see also ROIPOLY for creation of binary mask
%

out   = zeros(size(I,1), size(I,2)); % init error surface with zeros
J     = im2double(J);                % binary image mask

% if there exists no mask support, out is filled with 0's
if (sum(J(:)) <= 0), 
    return; 
end

fftJ  = fft2(J);        % fft of J
fftI  = fft2(I);        % fft of I (image mask)
fftIs = fft2(I.^2);     % fft of I sqaured
 
% symbolic constants for color channel weighting
RED_WEIGHT = 0.299;
GREEN_WEIGHT = 0.587;
BLUE_WEIGHT = 0.114;

% channel weights vector
w = [RED_WEIGHT, GREEN_WEIGHT, BLUE_WEIGHT];

% equation (2) from hts paper, augmented with
% (source) weighting (S) and normalization (hts thesis)
for color=1:3,
    out = out + w(color)*((real(ifft2(fftJ.*conj(data.fftSTs(:,:,color)))) - ...
        2 * real(ifft2(fftI(:,:,color).*conj(data.fftST(:,:,color)))) + ...
        real(ifft2(fftIs(:,:,color).*conj(data.fftS(:,:,color))))) ...
        ./ real(ifft2(fftJ.*conj(data.fftS(:,:,color)))));
end

% output error image must be rotated by 180 degrees (equiv. to mirror about both main axes)
% this ensures that each coordinate (x,y) in the error image stores the value for picking 
% the patch from the input texture with (x,y) as upper left bounding box cordinate.
% this is mainly to make the output more intuitive, as if we were circularly shifting
% the mask I and binary support J, where in reality we are circularly shifting the
% input texture T (see equation (1) in 'hybrid texture synthesis' paper)
out = rot90(out,2);
