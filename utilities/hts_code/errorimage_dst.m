function [out] = errorimage_dst(data,I,J)
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
% File errorimage_dst.m
%
%   [out] = errorimage_dst(data,I,J)
%
%   This subroutine will compute the error
%   surface for a given texure (stored in 'data'), image mask I and mask 
%   support function J. 
%
%   extended to DESTINATION IMPORTANCE WEIGHTING using a 
%   destination texture importance map D, described in hts paper/thesis:
%
%   'Hybrid Texture Synthesis', Nealen A.
%   Diplomarbeit (MSc thesis), Technische Universitaet Darmstadt
%
%   complexity O(N log N) with N = w x h
%
%   INPUT:
%   data - the input datastructure as defined in gen_input.m
%   I - the image mask
%   J - the binary mask support function
%
%   OUTPUT:
%   out - the error image with error values in [0,1]
%
%   The routine will return the error surface. indices are 1-based, meaning 
%   the value in the upper left corner is the error for shifting the mask 
%   by +1 pixel in x and y directions. the zero value for no translation is 
%   stored at index [max_height,max_width] (bottom right), which is eqivalent 
%   to no translation as the input texture is handled toroidally
%
%   see also ROIPOLY for creation of binary mask (in hybridsynthesizerec.m)
%

out = zeros(size(I,1), size(I,2)); % init error surface with zeros
J   = im2double(J);                % binary image mask

% if there exists no mask support, out is only 0's
if (sum(J(:)) <= 0), 
    return; 
end

% implementation of destination weighting by generating frequency map from I 
% (the destination overlap region). this must be done in each step as opposed 
% to the source frequency map which is only computed once on initialization.
D = frequencymap(I);

% for testing purposes. uniform weighting, so we can compare (debug) the metrics
% D = ones(size(I,1), size(I,2));

% convert to rgb
D = gray2rgb(D);
I = gray2rgb(I);
J = gray2rgb(J);

% cut off regions outside of J's support
JD = D .* J;

% other combinations we need for error metric
IJD  = I .* JD;
Is   = I.^2;
IsJD = Is .* JD;

% fourier transforms
fftJD  = fft2(JD);     % fft of JD
fftIJD = fft2(IJD);    % fft of IJD
 
% symbolic constants for color channel weighting
RED_WEIGHT = 0.299;
GREEN_WEIGHT = 0.587;
BLUE_WEIGHT = 0.114;

% channel weights vector
w = [RED_WEIGHT, GREEN_WEIGHT, BLUE_WEIGHT];

% equation (2) from hts paper, augmented with
% (destination) weighting (D) and normalization (hts thesis)
for color=1:3,
    out = out + w(color)*((sum(sum(IsJD(:,:,color))) - ...
        2 * real(ifft2(fftIJD(:,:,color).*conj(data.fftT(:,:,color)))) + ...
        real(ifft2(fftJD(:,:,color).*conj(data.fftTs(:,:,color))))) ...
        ./ sum(sum(JD(:,:,color))));
end

% output error image must be rotated by 180 degrees (equiv. to mirror about both main axes)
% this ensures that each coordinate (x,y) in the error image stores the value for picking 
% the patch from the input texture with (x,y) as upper left bounding box cordinate.
% this is mainly to make the output more intuitive, as if we were circularly shifting
% the mask I and binary support J, where in reality we are circularly shifting the
% input texture T (see equation (1) in 'hybrid texture synthesis' paper)
out = rot90(out,2);
