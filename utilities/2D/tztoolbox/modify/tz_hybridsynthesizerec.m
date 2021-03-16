function error = hybridsynthesizerec(data,overlap,patchlist)
%TZ_HYBRIDSYTHESIZEREC Modified from HYBRIDSYTHESIZEREC
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
% HYBRIDSYTHESIZEREC adaptive patch quilting with 'overlap resynthesis' 
%                    patch lapping improvement strategy. this is the recursive 
%                    part, thus the use of suffix REC.
%                    this is also the heart of the implementation
%
%   error = hybridsynthesizerec(data,overlap,patchlist)
%
%   This subroutine will synthesize a texture from an input
%   sample data.T, and use patches defined in patchlist.F and patchlist.V 
%   (faces and vertices) in the order they occur in the facelist F.
%   it extends image quilting to using a technique which we title
%   'image completion', this is implemented to reduce aliasing artifacts 
%   due to insufficient neighborhood matching along patch edges. also, this 
%   algorithm implements 'adaptive patch sampling', where we split a patch into
%   four patches when the patch-lapping error is above data.patch_delta_max
%
%   INPUT:
%   data      - the input 'data' as described in gen_input.m
%   overlap   - the number of overlap pixels
%   patchlist - an indexed face set structure with vertices V, faces F
%               and bounding box data BB (see gentoroidalquadmeshoffset.m
%               and genquadmeshoffset.m)
%
%   OUTPUT:
%   mse - the lapping-error for the patch encountered in this call
%

% global figure handle
global pos1;

% global result (out) and intermediate (pre) image
global out;
global pre;

% global source map (for ashikhmin/k-coherence search)
global source_map;

% k-nearest neighbors for each pixel in T (flattened, scanline-order)
% see hybridsynthesize.m and especially build_knn.m. here, we only need it to
% decide whether to use overlap_resynthesis_knn.m or overlap_resynthesis_exh.m
global knn;

% some global magic (pixelvalue) numbers, defined in hybridsynthesize.m
global NOT_YET_SYNTHESIZED;
global PER_FRAGMENT_SYNTH;

% the error, accumulated for each patch (in this call)
error_accum = 0;

% -----------------------------------------------------------------------------------------
% ALGORITHM STEP: iterate through all patches
% -----------------------------------------------------------------------------------------
% iterate through all output patches in order given by face list F
% note: the ordering has major impact on the synthesis results! we need
% causal neighborhoods as much as possible.
for posval=1:length(patchlist.F),

    % pos = patch number, for shifting start patch
    pos = posval;
   
    % -------------------------------------------------------------------------------------
    % ALGORITHM STEP: compute 'shift'
    % -------------------------------------------------------------------------------------
    % shift patch vertices (tmp data stored in 'patch') so that BB min of 
    % F(pos) is in position(1,1), toroidally
    % BB(pos,1) and BB(pos,2) are upper left BB corner row and col index
    % values in 'shift' will be <= 0 as the BB coords are positive integers >= 1
    % 'original_patch' is a list of patch vertices in clockwise ordering
    original_patch = [patchlist.V(patchlist.F(pos,:),1) patchlist.V(patchlist.F(pos,:),2)];
    shift = [1 1] - [patchlist.BB(pos,1) patchlist.BB(pos,2)];
    % -------------------------------------------------------------------------------------
    % ALGORITHM STEP: 'shift' patch so min bounding box coord is in (1 1)
    % -------------------------------------------------------------------------------------
    patch = circshiftpatch(original_patch,data.rows,data.cols,shift);

    % -------------------------------------------------------------------------------------
    % ALGORITHM STEP: compute binary support for ungrown patch
    % -------------------------------------------------------------------------------------
    % compute binary patch support PS (of size 'out') for the ungrown patch
    % note -> on 'patch-1': subtract 1 from patch coords (making them 0-based) to 
    % compensate for matlab's ROIPOLY fill conventions
    trnpatch = transpose(patch-1);
    PS = roipoly(out, [trnpatch(2,:)], [trnpatch(1,:)]);
    
    % -------------------------------------------------------------------------------------
    % ALGORITHM STEP: now shift patch and grow by dilation using user-defined 'overlap'
    %                 GPS stands for 'grown patch support' (binary)
    % -------------------------------------------------------------------------------------
    psize = max(patch) - 1;
    deltashift = [overlap overlap];
    GPS = circshift(PS, deltashift);
    GPS = imdilate(GPS, strel('square', 2*overlap+1));
    grownpatchsize = psize + 2*overlap;
    
    % some important error output. if the grown patch is larger than the 
    % input texture T, we must abort...
    if (grownpatchsize(1) > data.rows_t | grownpatchsize(2) > data.cols_t),
        disp('ERROR: grown patch is larger than input texture. aborting...');
        return;
    end
    
    % -------------------------------------------------------------------------------------
    % ALGORITHM STEP: 'deltashift' ungrown patch support 'PS' as well
    % -------------------------------------------------------------------------------------
    % also, shift PS by deltashift (so the ungrown patches in both PS and GPS are congruent)
    PS = circshift(PS,deltashift);
    
    % -------------------------------------------------------------------------------------
    % ALGORITHM STEP: shift copy of 'pre' and build binary support J and image mask I
    % -------------------------------------------------------------------------------------
    % circularly shift temporary copy of (globally defined) 'pre' by (shift + deltashift)
    % note: pre is partially synthesized result with values of 
    % (NOT_YET_SYNTHESIZED 0 0) in pixels that are yet to be synthesized
    pre_shift = circshift(pre, shift+deltashift);
    [I,J] = buildimagemask(size(data.T), grownpatchsize, GPS, pre_shift);
    
    % -------------------------------------------------------------------------------------
    % ALGORITHM STEP: compute error image from I,J, fftS, fftST and fftSTs
    % -------------------------------------------------------------------------------------
    % now we have I, J and all other data (in 'data') and can generate error image for 
    % the image mask to locate the best match (remember, reference point is upper left 
    % corner of support function area and image mask bounding box)
        
    % switch on 'metric' string in 'data' (see gen_input.m)
    switch data.metric
    case 'src'
        errorimg = errorimage_src(data, I, J);
    case 'dst'
        errorimg = errorimage_dst(data, I, J);
    case 'sum'
        errorimg = errorimage_sum(data, I, J);
    case 'sqdiff'
        errorimg = errorimage_sqdiff(data, I, J);
    otherwise
        % this is equivalent to setting metric to 'simple'
        errorimg = errorimage(data, I, J);
    end

    % NOTE: value in max(errorimage) (bottom right corner) is error for picking patch with 
    % BB_grown_min at (1,1) (no shift)

    % -------------------------------------------------------------------------------------
    % ALGORITHM STEP: trim errorimg if input T doesn't wrap 
    % -------------------------------------------------------------------------------------
    % trim bottom and right strip of err (if T doesn't wrap) -> (NOTE: 'wrap' CAN 
    % cause excessive verbatim copying if we do not introduce countermeasures)
    if (strcmp('nowrap',data.mode)),
        errorimg = errorimg(1:size(errorimg,1)-grownpatchsize(1)+1,...
            1:size(errorimg,2)-grownpatchsize(2)+1);
    end
        
    % -------------------------------------------------------------------------------------
    % ALGORITHM STEP: find some candidates for picking with low error, but diverse
    % -------------------------------------------------------------------------------------
    % find some candidates depending on user defined 'errtol' and choose one at random
    % must use abs() because rounding errors (can) lead to negative, close-to-zero numbers. 
    min_error = min(errorimg(:));
    minimumerror = abs(min_error);
    error_tolerance = (1 - minimumerror) * data.errtol;
    [jbest, ibest] = find(errorimg <= (error_tolerance+minimumerror));
    
    % choose one at randindex
    randindex = ceil(rand * length(jbest));
    % bestpos is the upper left bounding box corner of best GROWN patch
    bestpos = [jbest(randindex) ibest(randindex)];
    
    % uncomment the following 4 lines to hardcode placement of the first (seed)
    % patch. setting data.errtol to 0 and hardcoding the first patch yields
    % a deterministic, reproduceable synthesis result

    % if (min_error == 0),
    %  disp('warning: hardcoded first patch!!!');
    %  bestpos = [128 128];        
    % end
    
    % -------------------------------------------------------------------------------------
    % ALGORITHM STEP: shift T by bestpos and compute the error surface S between I and T
    %                 mse -> mean squared error, ssd -> sum of square differences
    % -------------------------------------------------------------------------------------
    shiftT = circshift(data.T,-bestpos);
    [S,bestmatch_mse,bestmatch_ssd,numpixels] =...
        errorsurface(shiftT,grownpatchsize(1),grownpatchsize(2),J,I);

    % use mean squared error (or alternatively 'metric' error) to determine 
    % when to split patches
    patch_overlap_error = bestmatch_mse;
    % patch_overlap_error = min_error;
    
    % we currently stop splitting at patchsize defined in 'data' input structure
    finalsplit = 2*data.min_patchsize;
    if (patch_overlap_error > data.patch_delta_max & ...
            psize(1) >= finalsplit & psize(2) >= finalsplit),
        % ---------------------------------------------------------------------------------
        % ALGORITHM STEP: if patch_overlap_error too large, 
        %                 split (if possible), and recurse
        % ---------------------------------------------------------------------------------
        % make four patches out of this one patch by creating five new vertices and
        % four faces from the 9 vertices, see GenUniformQuadMeshWithOffset.m
        % NOTE: this implementation only works with regular quad patches!!!
        row_offset = patchlist.BB(pos,1) - 1;
        col_offset = patchlist.BB(pos,2) - 1;
        row_psize = psize(1)/2;
        col_psize = psize(2)/2;
        % also split pixel overlap for the next level of patches. NOTE:
        % we want a minimum of 3 lapping pixels (equivalent to a 7x7 nbhd for 1 pixel)
        overlap_new = max(3,ceil(overlap/2));

        [splitpatchlist.V,splitpatchlist.F,splitpatchlist.BB] = ...
            genquadmeshoffset(row_psize,col_psize,2,2,row_offset,col_offset);
        sub_error = hybridsynthesizerec(data,overlap_new,splitpatchlist);
        
        % accumulate the mse for this call
        error_accum = error_accum + sub_error;

        % some output
        if (data.verbatim),
            disp(sprintf('splitting. error: %f, patchsize: %d->%d, sub_error: %f',...
                bestmatch_mse, row_psize*2, row_psize, sub_error));
        end
    else
        
        % accumulate the error for this call
        error_accum = error_accum + bestmatch_mse;
        
        % ---------------------------------------------------------------------------------
        % ALGORITHM STEP: build alpha blending mask to feather valid pixels
        % ---------------------------------------------------------------------------------
        numalphasteps = overlap; % blend over entire overlap width
        if (data.feather == 0), numalphasteps = 0; end
        alphamask = buildalphamask(PS,J,grownpatchsize,numalphasteps);

        % ---------------------------------------------------------------------------------
        % ALGORITHM STEP: build validpixel mask and 'J OR PS' surface (=PSall)
        %                 PSall is original patch + area of grown patch which overlaps
        %                 already synthesized pixels
        % ---------------------------------------------------------------------------------
        % now we must use the entire patch plus the mask area as the lapping mask, set 
        % the validpixel-mask to 1 in errorvalid regions, and blackout invalid regions 
        % for later filling. note: a setting of 0 for data.pixel_delta_max results in 
        % simple feathering (if data.feather == 1), or abrupt borders
        PSall = PS;
        validpixels = PS;
        for jj=1:grownpatchsize(1),
            for ii=1:grownpatchsize(2),
                if (J(jj,ii) == 1 & PS(jj,ii) == 0),
                    PSall(jj,ii) = 1;
                    if (S(jj,ii) <= data.pixel_delta_max), validpixels(jj,ii) = 1; else
                        validpixels(jj,ii) = 0; end
                end
            end
        end
        
        % ---------------------------------------------------------------------------------
        % ALGORITHM STEP: build traversal map for the per-fragment resynthesis in J
        % ---------------------------------------------------------------------------------
        pixeltraversalmap = buildtraversalmap(GPS,PS,J,validpixels,grownpatchsize);
    
        % ---------------------------------------------------------------------------------
        % ALGORITHM STEP: compose the best match with the output and pre images
        %                 and prepare 'out' and 'pre' for per-fragment synthesis
        % ---------------------------------------------------------------------------------
        % use empty image ('bestpatch') of same size as output image, paste the best 
        % patch from texture example in upper left, circularly shift it by values stored in 
        % 'shift+deltashift' and add it to 'out' image (and 'pre' image)
    
        % build a copy of (flattened) indices into T for the source_map by
        % shifting a copy of the index_map so it's congruent with shiftT
        shift_index_map = circshift(data.index_map,-bestpos);
        best_src_indices = zeros(data.rows, data.cols);
        
        % first clear the 'bestpatch' image again
        bestpatch = zeros(data.rows, data.cols, 3); 
        % build best patch in 'bestpatch'. concurrently build the source indices.
        for jj=1:grownpatchsize(1),
            for ii=1:grownpatchsize(2),
                % if we have patch support here...
                if (PSall(jj,ii) == 1), 
                    % ...add pixel from shifted T to 'bestpatch'. 
                    % note: shiftT was already used to compute error surface S
                    bestpatch(jj,ii,1:3) = shiftT(jj,ii,1:3); 
                    % and add source index
                    best_src_indices(jj,ii) = shift_index_map(jj,ii);
                end
            end
        end
    
        % shift this best patch 'bestpatch' to final destination and compose with 
        % 'out' and 'pre' based on pixel validity and alpha values
        bestpatch = circshift(bestpatch,-(shift+deltashift));
        best_src_indices = circshift(best_src_indices,-(shift+deltashift));
        PSall = circshift(PSall,-(shift+deltashift));
        validpixels = circshift(validpixels,-(shift+deltashift));
        alphamask = circshift(alphamask,-(shift+deltashift));
        pixeltraversalmap = circshift(pixeltraversalmap,-(shift+deltashift));
  
        % for all pixels within support of PS or J (= PSall)
        for jj=1:data.rows,
            for ii=1:data.cols,
                % replace pixel values where validpixels is 1, set to value 
                % (PER_FRAGMENT_SYNTH 0 0) where validpixels is 0 (in pre)
                % reminder: PSall is ungrown patch PLUS mask support (overlap)
                if (PSall(jj,ii) == 1),
                    ispixelvalid = validpixels(jj,ii);
                    alpha = alphamask(jj,ii);
                    if (ispixelvalid == 1),
                        % VALID PIXEL
                        % compose final pixel into out and pre, using basic OVER compositing
                        out(jj,ii,1:3) = (1-alpha)*out(jj,ii,1:3) + alpha*bestpatch(jj,ii,1:3);
                        pre(jj,ii,1:3) = (1-alpha)*pre(jj,ii,1:3) + alpha*bestpatch(jj,ii,1:3); 
                        
                        % build source_map. if alpha > 0.5, use source_index from 
                        % new patch, otherwise retain old source value (= do nothing)
                        if (alpha > 0.5),
                            source_map(jj,ii) = best_src_indices(jj,ii);
                        end
                    
                    else
                        % INVALID PIXEL
                        % initialize out and pre for per-fragment synthesis
                        % pre now has (NOT_YET_SYNTHESIZED 0 0) where no synthesis has 
                        % taken place, and (PER_FRAGMENT_SYNTH 0 0) where pixels are to 
                        % be resynthesized per-fragment in this step 
                        out(jj,ii,1:3) = 0; 
                        % use (red) color gradient for displayed, intermediate result
                        factor = 1/pixeltraversalmap(jj,ii);
                        out(jj,ii,1) = min(factor*2.5,1);            
                        pre(jj,ii,2:3) = 0; pre(jj,ii,1) = PER_FRAGMENT_SYNTH;
                    end
                end
            end
        end
        
        % ---------------------------------------------------------------------------------
        % ALGORITHM STEP: fill invalid pixels using the traversal map
        % ---------------------------------------------------------------------------------
        % fill holes where invalid lapping pixels exist in traversal order 

        % if we pass in a knn-datastructure (globally defined, see hybridsynthesize.m)
        % for T, use it. otherwise default to exhaustive (fourier-based) search.
        if (knn ~= 0),
            % do ashikhmin/k-coherence search
            % using box-shaped neigborhood defined in 'data'
            overlap_resynthesis_knn(data,pixeltraversalmap);
        else
            % perform exhaustive, fft-based search 
            % using box-shaped neigborhood defined in 'data'
            overlap_resynthesis_exh(data,pixeltraversalmap);
        end

        % always display the current result (patch-progress)
        figure(1); imshow(out,'truesize')
    end
    
    % mse for this call: divide accumulated mse by the number of patches in this call
    error = error_accum/length(patchlist.F);
end
