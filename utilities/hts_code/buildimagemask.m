function [I,J] = buildimagemask(size,grownpatchsize,support,prevsynth);
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
% File buildimagemask.m
%   This subroutine returns the image mask I and binary support of size 'size'
%
%   [I,J] = buildimagemask(size,grownpatchsize,support,prevsynth)
%
%   INPUT:
%   size - [r c] size of I and J
%   grownpatchsize - [r c] size of the grown patch for iteration
%   support - the 2D binary support function of the grown patch (GPS)
%   prevsynth - the 2D already synthesized parts of the result, with
%               (NOT_YET_SYNTHESIZED 0 0) where no synthesis has yet taken place
%
%   OUTPUT:
%   I - the image mask (rgb values)
%   J - binary support for I (needed in error computation)
%

% some global magic (pixelvalue) numbers
global NOT_YET_SYNTHESIZED;

% check the rectangular area of the grown patch (BBmin_gp, BBmax_gp) for values 
% where prevsynth is not equal to init value (NOT_YET_SYNTHESIZED 0 0) and support is 1. 
% build image mask (I) and mask support (J) from these values (I and J must have the 
% size of the input texture for error image computation)
I = zeros(size(1), size(2), 3);
J = zeros(size(1), size(2));
for jj=1:grownpatchsize(1),
    for ii=1:grownpatchsize(2),
        if ((support(jj,ii) == 1) & (prevsynth(jj,ii,1) ~= NOT_YET_SYNTHESIZED)), 
            % if the grown patch has support here, and there exists a valid pixel value
            % add support for this pixel (J)
            J(jj,ii) = 1;
            % and store the pixel value in the image mask (I)
            I(jj,ii,1:3) = prevsynth(jj,ii,1:3);
        end
    end
end
