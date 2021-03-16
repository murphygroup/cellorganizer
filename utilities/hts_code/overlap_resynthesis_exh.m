function overlap_resynthesis_exh(data,traversalmap)
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
% File overlap_resynthesis_exh.m
%   This subroutine will fill holes in an input texture
%   on a per-fragment basis by exhaustively searching in 'data.T'
%
%   overlap_resynthesis_exh(data,traversalmap)
%
%   INPUT:
%   data         - the input datastructure as defined in gen_input.m
%   traversalmap - the traversalorder for pixel synthesis (integer values 
%                  starting at 2)
%
%   OUTPUT:
%   this function modifies globally defined image arrays, see hybridsynthesize.m
%

% global figure handle
global pos1;

% global result and intermediate images
global out;
global pre;

% some global magic (pixelvalue) numbers, defined in hybridsynthesize.m
global NOT_YET_SYNTHESIZED;
global PER_FRAGMENT_SYNTH;

% abort, if traversalmap is made of only 1's (= all pixels synthesized)
if (max(traversalmap(:)) < 2), return; end

overlap = data.exh_overlap;

% for each traversal level in the traversalmap
numtraversallevels = max(traversalmap(:));
for traversallevel=2:numtraversallevels,
    
    % optional: display per-level progress
    if (data.display_resynthesis),
        figure(1); imshow(out)
    end
    
    % iterate through all output pixels in scanline order
    for j=1:size(out,1),     % j = rows
        for i=1:size(out,2), % i = columns
            
            % only perform synthesis for invalid pixels in current traversallevel
            if (traversalmap(j,i) == traversallevel),
                
                % xpos, ypos: coordinate (1-based) of currently synthesized pixel
                ypos = j; % row index
                xpos = i; % col index
                
                % compute upper left corner of patch, grown by overlap (ymask, xmask) and
                % shift the 'pre' image so that this 'point of reference' is in upper left 
                % corner (shift_pre), index(1,1) (this guarantees wrapping textures)
                yshiftmask = ypos-overlap-1; xshiftmask = xpos-overlap-1;
                shift_pre = circshift(pre, [-yshiftmask -xshiftmask]);
                
                % now check the area of the grown pixel for 'valid' values in 'pre'
                % and build image mask (I) and mask support (J) from them (I and J must
                % have the size of the input texture for error computation)
                xmax = 2*overlap + 1; ymax = 2*overlap + 1;
                I = zeros(data.rows_t, data.cols_t, 3);
                J = zeros(data.rows_t, data.cols_t); % binary support for image mask
                for jj=1:ymax,
                    for ii=1:xmax,
                        % if there exists a valid pixel value (implicitly in the red channel
                        % see hybridsynthesize.m for more info on this)
                        if (shift_pre(jj,ii,1) ~= NOT_YET_SYNTHESIZED &...
                            shift_pre(jj,ii,1) ~= PER_FRAGMENT_SYNTH),
                            % add support for this pixel
                            J(jj,ii) = 1;
                            % and store the pixel value in the image mask (I) 
                            I(jj,ii,1:3) = shift_pre(jj,ii,1:3);
                        end
                    end
                end
                % now we have I, J and T and can generate error surface for the image mask to
                % locate the best match (remember, reference point is upper left corner of 
                % support function area and image mask
                err = errorimage(data, I, J);
                
                % shift error image 'err' so that upper left corner is error value 
                % for picking pixel at that position
                err = circshift(err,[overlap+1 overlap+1]);
                % trim lower and right edge, implying that T doesn't wrap
                err = err(1:size(err,1)-overlap, 1:size(err,2)-overlap);
                % find some candidates and choose one at random
                % must use abs() because rounding errors lead to negative
                % close-to-zero numbers. pick only the best match for each pixel
                [jbest, ibest] = find(err <= 1.00*abs(min(err(:))));
                c = ceil(rand * length(jbest));
                
                % bestypos, bestxpos - pixel coords of best match in T
                bestypos = jbest(c); bestxpos = ibest(c);

                % paste that pixel into 'out' and 'pre'
                pixel = data.T(bestypos, bestxpos, 1:3);
                out(ypos, xpos, 1:3) = pixel;
                pre(ypos, xpos, 1:3) = pixel;
                
                % then blend the four neighbors with the existing synthesis result
                % do not let the texture wrap in this step
                blendpix = [];
                if (bestypos < data.rows_t), blendpix = cat(3, blendpix, [+1 0]); end
                if (bestypos > 1),           blendpix = cat(3, blendpix, [-1 0]); end
                if (bestxpos < data.cols_t), blendpix = cat(3, blendpix, [0 +1]); end
                if (bestxpos > 1),           blendpix = cat(3, blendpix, [0 -1]); end
                
                for bpixel=1:size(blendpix,3),
                    % xn, yn - pixel coords of current pixel neighbors in out and pre
                    yn = ypos + blendpix(1,1,bpixel); xn = xpos + blendpix(1,2,bpixel);
                    if (yn > 0 & xn > 0 & yn < size(out,1) & xn < size(out,2)),
                        % if synthesis result already exists for this neighbor...
                        if (pre(yn, xn, 1) ~= NOT_YET_SYNTHESIZED &...
                            pre(yn, xn, 1) ~= PER_FRAGMENT_SYNTH),
                            % ...blend into out and pre with alpha of 0.5
                            % ynT, xnT - pixel coords of chosen pixel neigbors in T
                            ynT = bestypos + blendpix(1,1,bpixel);
                            xnT = bestxpos + blendpix(1,2,bpixel);
                            out(yn, xn, 1:3) = 0.5 * out(yn, xn, 1:3) + ...
                                0.5 * data.T(ynT, xnT, 1:3);
                            pre(yn, xn, 1:3) = 0.5 * out(yn, xn, 1:3) + ... 
                                0.5 * data.T(ynT, xnT, 1:3);
                        end                    
                    end
                end
            end
        end
    end
end
