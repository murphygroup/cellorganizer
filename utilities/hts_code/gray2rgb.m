function [out] = gray2rgb(I)
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
% GRAY2RGB simply transfer grayscale image into rgb mode, if
%          the image is not already rgb
%
%   [out] = gray2rgb(I)
%

out = I;
if (length(size(I)) > 2), return; end 

out = zeros(size(I,1), size(I,2), 3);
for j=1:size(I,1),
    for i=1:size(I,2),
        out(j,i,1:3) = I(j,i);
    end
end
