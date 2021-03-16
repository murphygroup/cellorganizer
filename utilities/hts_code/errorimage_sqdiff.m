function [out] = errorimage_sqdiff(data,I,J)
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
% File errorimage_sqdiff.m
%
%   [out] = errorimage_sqdiff(data,I,J)
%
%   This subroutine will compute the error
%   surface for a given texure (stored in 'data'), image mask I and mask 
%   support function J. 
%
%   extended to bilateral (source/destination)
%   importance weighting from hybrid texture synthesis thesis
%
%   'Hybrid Texture Synthesis', Nealen A.
%   Diplomarbeit (MSc thesis), Technische Universitaet Darmstadt
%
%   using sqdiff(S,D) as the importance map compositing function
%
%   complexity O(N log N) with N = w x h
%
%   INPUT:
%   data - the input datastructure as defined in hybridsynthesize.m and gen_input.m
%   I    - the image mask
%   J    - the binary mask support function
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

out = zeros(size(I,1), size(I,2)); % init error surface with zeros
J   = im2double(J);                % binary image mask

% if there exists no mask support, out is filled with 0's
if (sum(J(:)) <= 0), 
    return; 
end

fftJ  = fft2(J);        % fft of J
fftI  = fft2(I);        % fft of I (image mask)
Is = I.^2;
fftIs = fft2(Is);       % fft of I sqaured

% implementation of destination weighting by generating frequency map from I 
% (the destination overlap region). this must be done in each step as opposed 
% to the source frequency map which is only computed once on initialization.
D = frequencymap(I);

% for testing purposes. uniform weighting, so we can compare (debug) the metrics
% D = ones(size(I,1), size(I,2));

% add a little epsilon to get away from excessive zeros in sqdiff(S,D). 
% especially when both importance maps are set to all 1's (uniform weighting)
D = D + 0.001;

% cut off regions outside of J's support
D = J .* D;
% convert D to rgb
D = gray2rgb(D);

% things we need for error computation (related to destination frequency)
Ds = D.^2;
DsI = Ds.*I;
DsIs = Ds.*Is;
DIs = D.*Is;
DI = D.*I;
fftD = fft2(D);
fftDI = fft2(DI);
fftDs = fft2(Ds);
fftDsI = fft2(DsI);
fftDIs = fft2(DIs);

% symbolic constants for color channel weighting
RED_WEIGHT = 0.299;
GREEN_WEIGHT = 0.587;
BLUE_WEIGHT = 0.114;

% channel weights vector
w = [RED_WEIGHT, GREEN_WEIGHT, BLUE_WEIGHT];

% equation (2) from hts paper, augmented with 
% source/destination weighting and normalization (from 'hybrid texture synthesis' 
% thesis)
for color=1:3,
    out = out + ...
        w(color)*(...
        (real(ifft2(fftIs(:,:,color).*conj(data.fftSs(:,:,color)))) - ...
        2 * real(ifft2(fftI(:,:,color).*conj(data.fftTSs(:,:,color)))) + ...
        real(ifft2(fftJ.*conj(data.fftTsSs(:,:,color)))) - ...
        2 * real(ifft2(fftDIs(:,:,color).*conj(data.fftS(:,:,color)))) + ...
        4 * real(ifft2(fftDI(:,:,color).*conj(data.fftST(:,:,color)))) - ...
        2 * real(ifft2(fftD(:,:,color).*conj(data.fftSTs(:,:,color)))) - ...
        2 * real(ifft2(fftDsI(:,:,color).*conj(data.fftT(:,:,color)))) + ...
        real(ifft2(fftDs(:,:,color).*conj(data.fftTs(:,:,color)))) + ...
        sum(sum(DsIs(:,:,color)))) ...
        ./ (real(ifft2(fftJ.*conj(data.fftSs(:,:,color)))) - ...
        2*real(ifft2(fftD(:,:,color).*conj(data.fftS(:,:,color)))) + sum(sum(Ds(:,:,color))))...
        );
end

% output error image must be rotated by 180 degrees (equiv. to mirror about both main axes)
% this ensures that each coordinate (x,y) in the error image stores the value for picking 
% the patch from the input texture with (x,y) as upper left bounding box cordinate.
% this is mainly to make the output more intuitive, as if we were circularly shifting
% the mask I and binary support J, where in reality we are circularly shifting the
% input texture T (see equation (1) in 'hybrid texture synthesis' paper)
out = rot90(out,2);
