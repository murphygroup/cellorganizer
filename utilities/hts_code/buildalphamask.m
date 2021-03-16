function [alphamask] = buildalphamask(PS,J,grownpatchsize,steps);
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
% File buildalphamask.m
%   This subroutine returns the alphamask for the input patch 
%   with binary support PS
%
%   [alphamask] = buildmask(PS,J,grownpatchsize,steps)
%
%   INPUT:
%   PS - 2D binary support function for the ungrown patch
%   J  - 2D image mask (overlap region) binary support function
%   grownpatchsize - [r c] size of the grown patch for iteration
%                    a two-element vector
%   steps - number of alpha blending steps
%
%   OUTPUT:
%   alphamask - the alpha blending mask for this patch
%

SE = strel('square',3); 
PSdilated = PS; PSdilatedprev = PS;
alphamask = im2double(PS);
% 'steps' defines the alpha-levels
if (steps > 0), alphavalues = 0:1/steps:1-(1/steps); end;
for dilationstep=1:steps,
    PSdilated = imdilate(PSdilated,SE);
    % set 1's to 0 where there exists no image mask support 
    % (J == 0 && PS == 0, only allow dilation into mask)
    for jj=1:grownpatchsize(1),
        for ii=1:grownpatchsize(2),
            if ((J(jj,ii) == 0 & PS(jj,ii) == 0)), PSdilated(jj,ii) = 0; end
        end
    end
    % now get the difference to previous step, set alpha mask, and prepare next iteration
    diff = (PSdilated & ~PSdilatedprev);
    for jj=1:grownpatchsize(1),
        for ii=1:grownpatchsize(2),
            if (diff(jj,ii) == 1), alphamask(jj,ii) = alphavalues(steps-dilationstep+1); end
        end
    end
    PSdilatedprev = PSdilated;
end
