function traversalmap = buildtraversalmap(GPS,PS,J,validpixels,grownpatchsize);
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
% File buildtraversalmap.m
%   This subroutine returns the traversalmap for overlap-resynthesis
%
%   traversalmap = buildtraversalmap(GPS,PS,J,validpixels,grownpatchsize)
%
%   INPUT:
%   GPS - 2D binary support function for the grown patch
%   PS  - 2D binary support function for the ungrown patch
%   J   - 2D image mask (lapping area) binary support function
%   validpixels - a 2D array with valid pixels set to 1, otherwise 0
%   grownpatchsize - [r c] size of the grown patch for iteration
%
%   OUTPUT:
%   traversalmap - the traversal levels for per-fragment synthesis with 1's 
%                  as already synthesized pixels, 2's are the next level
%                  3's the next... etc. max(traversalmap(:)) returns the
%                  highest level in the map
%

% we define a traversal order for the invalid pixel overlap-resynthesis.
% a straightforward method is to iteratively dilate the valid pixels with a small
% structuring element and set a key value in each step. this value-array must be passed
% to the pixel resynthesizer (overlap_resynthesis). the dilation occurs until the 
% entire patch (PS OR J) is filled and the traversal order is set in 'traversalmap', 
% where 1's are already synthesized pixels, 2's are the next level, the 3's, etc...

% note: shift circularly by sm pixels for dilation. shift back afterwards
sm = [0 0];
% optionally: dilate from new patch AND existing synthesis result
%             note that in this case, sm must be set to [1 1]
% traversalmap = im2double(validpixels | ~GPS);
traversalmap = im2double(validpixels); % only traverse from new patch + validpix
% uncomment the following line is sm is not [0 0]
% traversalmap = circshift(traversalmap, sm);
dilatevalidpixels = traversalmap;
dilatevalidpixelsprev = traversalmap;
SE = strel('disk',1);
% check if pixels must be inserted into the traversal map (openpixels > 0)
openpixels = 0;
for jj=1:grownpatchsize(1),
    for ii=1:grownpatchsize(2),
        if ((PS(jj,ii) == 1 | J(jj,ii) == 1) & dilatevalidpixels(jj+sm(1),ii+sm(2)) == 0),
            openpixels = openpixels + 1;
        end
    end
end
% start at level 2, level 1 represents already valid pixels
traversallevel = 2;
while (openpixels > 0),
    % dilate the current validpixels image and compute the difference to the previous
    dilatevalidpixels = imdilate(dilatevalidpixels,SE);
    traversaldiff = (dilatevalidpixels & ~dilatevalidpixelsprev);
    openpixels = 0;
    for jj=1:grownpatchsize(1),
        for ii=1:grownpatchsize(2),
            if ((PS(jj,ii) == 1 | J(jj,ii) == 1) & dilatevalidpixels(jj+sm(1),ii+sm(2)) == 0), 
                openpixels = openpixels + 1; 
            end
            % set the current level in the traversal map for newly reached pixels
            % within the image mask support J
            if ((PS(jj,ii) == 0 & J(jj,ii) == 1) & traversaldiff(jj+sm(1),ii+sm(2)) == 1), 
                traversalmap(jj+sm(1),ii+sm(2)) = traversallevel;
            end
        end
    end
    % increment the traversal level and prepare next step
    traversallevel = traversallevel + 1;
    dilatevalidpixelsprev = dilatevalidpixels;
end

% uncomment the following line is sm is not [0 0]
% traversalmap = circshift(traversalmap, -sm);
