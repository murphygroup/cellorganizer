function patch = circshiftpatch(P,rows,cols,shift)
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
% CIRCSHIFTPATCH circularly shift a patch based on pixel rows and cols
%
%   patch = circshiftpatch(P,rows,cols,shift)
%
%   simple circular shift of a patch with integer (pixel) vertex coordinates
%
%   input:
%   P  - the input patch, a list of 2D vertices
%   rows  - num pixel rows in final image
%   cols  - num pixel cols in final image
%   shift - [rows cols] array (1x2) with the shift values
%
%   output:
%   patch - the shifted patch
%
%   NOTE:
%   used in hybridsynthesizerec.m, only works for negative shift values
%   do not use outside of the context of this package!
%

patch = [(P(:,1) + shift(1)) (P(:,2) + shift(2))];
% number of rows in 'patch' defines numvertices of this patch
% note: size(patch,1) equals the number of patch vertices (=n for n-gon)
numverts = size(patch,1);
for n=1:numverts,
    % if anything has been shifted out of bounds, wrap it around
    if (patch(n,1) < 1), patch(n,1) = rows + patch(n,1); end % row index
    if (patch(n,2) < 1), patch(n,2) = cols + patch(n,2); end % col index
end
