function fmap = frequencymap(in)
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
% FREQUENCYMAP return a map in [0.2,1] where 1 represents high frequency
%              content, and 0.2 low frequency content
%
%   fmap = frequencymap(in)
%
%   input:
%   in - the input image
%
%   output:
%   fmap - the frequencymap in [0.2,1]
%
%   NOTE:
%   simply uses a high-pass filter. used to identify texture features.
%   see error metrics, which weight source and/or destination frequency-maps
%   using various compositing operators (errorimage_src.m, errorimage_dst.m,
%   errorimage_sum.m and errorimage_sqdiff.m)
%
%   gen_input.m also makes use of frequencymap.m to generate the source
%   importance map S
%

cutoff = 0.1;
order = 41;
padding = 10;
origsize = size(in);
startindex = padding + 1;
endindex_row = startindex + origsize(1) - 1;
endindex_col = startindex + origsize(2) - 1;

% convert to grayscale if necessary
if (length(size(in)) > 2),
    in = rgb2gray(in);
end
in = padarray(in,[padding padding],'symmetric','both');

% construct a 2D mesh of 'order' dimensions in the interval [0,1]
[f1,f2] = freqspace(order,'meshgrid');

% construct 2D tophat, high-pass filter in frequency space
Hd = zeros(order);
d = find(f1.^2+f2.^2 < cutoff^2);
% low-pass top-hat
Hd(d) = 1;
% make it a high-pass 'bucket'
Hd = Hd-1;

% sample from frequency space to produce 
% image space filter coefficients
h = fsamp2(Hd);

% and finally, apply this filter
in_f = filter2(h,in);

% map values to interval [0.2 1]

% push to 0 (no negative values)
in_f = in_f - min(in_f(:));

if (sum(in_f(:)) > 0),
    % map to [0 0.8]
    in_f = in_f ./ (1.25*max(in_f(:)));
    % push to [0.2 1]
    in_f = in_f + 0.2;
end

% cut away the padding
fmap = in_f(startindex:endindex_row, startindex:endindex_col, :);
