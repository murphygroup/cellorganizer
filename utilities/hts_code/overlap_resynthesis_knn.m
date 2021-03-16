function overlap_reynthesis_knn(data,traversalmap)
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
% File overlap_resynthesis_knn.m
%   This subroutine will fill holes in an input texture
%   on a per-fragment basis by ashikhmin/k-coherence search in 'data.T'
%
%   overlap_resynthesis_knn(data,traversalmap)
%
%   INPUT:
%   data         - the input datastructure as defined in gen_input.m
%   traversalmap - the traversalorder for pixel synthesis (integer values 
%                  starting at 2)
%   this function uses the globally defined knn datastructure and source_map, 
%   see hybridsynthesize.m
%
%   OUTPUT:
%   this function modifies globally defined image arrays, see hybridsynthesize.m
%

% global figure handle
global pos1;

% global result and intermediate images
global out;
global pre;

% global source map (for ashikhmin/k-coherence search)
global source_map;

% k-nearest neighbors for each pixel in T (flattened indices, scanline-order)
global knn;

% symbolic constants for color channel weighting
RED_WEIGHT = 0.299;
GREEN_WEIGHT = 0.587;
BLUE_WEIGHT = 0.114;

% channel weights vector
w = [RED_WEIGHT, GREEN_WEIGHT, BLUE_WEIGHT];

% some magic (pixelvalue) numbers
global NOT_YET_SYNTHESIZED;
global PER_FRAGMENT_SYNTH;

overlap = data.knn_overlap;

% abort, if traversalmap is made of only 1's (= all pixels synthesized)
if (max(traversalmap(:)) < 2), return; end

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

                % compute upper left corner of patch, grown by overlap 
                % (yshiftmask, xshiftmask) and shift the 'pre' image so that this 
                % 'point of reference' is in upper left corner, index(1,1) 
                % (this guarantees wrapping textures)
                row_shift = ypos-overlap-1;
                col_shift = xpos-overlap-1;
                shift_pre = circshift(pre, [-row_shift -col_shift]);
                shift_source_map = circshift(source_map, [-row_shift -col_shift]);
                
                % (initially empty) candiate list of pixels for L2 comparison 
                % with the target pixel
                candidate_list = [];
                
                % now check the area of the grown pixel for 'valid' values in 'preout'
                r_max = 2*overlap + 1;
                c_max = 2*overlap + 1;
                mid_coord = overlap + 1;
                for r=1:r_max,
                    for c=1:c_max,
                        % if there exists a valid pixel value ...
                        if (shift_pre(r,c,1) ~= NOT_YET_SYNTHESIZED &...
                            shift_pre(r,c,1) ~= PER_FRAGMENT_SYNTH),
                            % ...consider this pixels properly shifted pixel for 
                            % the candidate_list (see ashikhmin's paper on this)
                            source_index = shift_source_map(r,c);
                            index_map_r = floor((source_index-1)/data.cols_t)+1;
                            index_map_c = source_index - (index_map_r - 1)*data.cols_t;
                            % the values for shifting in index_map
                            shift_r = mid_coord - r;
                            shift_c = mid_coord - c;
                            % shifted indices for candidate in index_map
                            % IMPORTANT: make sure they wrap!!!
                            index_map_r = add_and_wrap(index_map_r,shift_r,data.rows_t);
                            index_map_c = add_and_wrap(index_map_c,shift_c,data.cols_t);
                            % current candidate_index...
                            candidate_index = data.index_map(index_map_r,index_map_c);
                            % ... used to index into k-nearest neighbors of candidate_index
                            if (data.k > size(knn,2)),
                                % if data.k is larger than available set, use all ...
                                new_candidates = ...
                                    cat(2,candidate_index,knn(candidate_index,:));
                            elseif (data.k > 1),
                                % ... otherwise limit to k-1 candidates from the knn set
                                %     plus the current candidate_index (=k candidates)
                                new_candidates = ...
                                    cat(2,candidate_index,knn(candidate_index,1:(data.k-1)));
                            else
                                % set k=1 for pure ashikhmin search
                                new_candidates = candidate_index;
                            end
                            
                            % now add these candidates to candidate_list
                            candidate_list = cat(2, candidate_list, new_candidates);
                        end
                    end
                end
                
                % make each candidate unique (remove duplicates)
                candidate_list = unique(candidate_list);
                
                % now compare each candidate in the list with the target neighborhood
                % and keep the index of the bestmatch
                lowest_error = 1000000000000;
                for L=1:length(candidate_list),
                    candidate_index = candidate_list(L);
                    row_t = floor((candidate_index-1)/data.cols_t)+1;
                    col_t = candidate_index - (row_t - 1)*data.cols_t;
                    row_shift = row_t-overlap-1; 
                    col_shift = col_t-overlap-1;
                    % compare neighborhoods in valid regions
                    current_error = 0;
                    for r=1:r_max,
                        row_index = add_and_wrap(row_shift,r,data.rows_t);
                        for c=1:c_max,
                            col_index = add_and_wrap(col_shift,c,data.cols_t);
                            % if there exists a valid pixel value ...
                            if (shift_pre(r,c,1) ~= NOT_YET_SYNTHESIZED &...
                                shift_pre(r,c,1) ~= PER_FRAGMENT_SYNTH),
                                % ... add this pixel to the error
                                for col=1:3,
                                    current_error = current_error + ...
                                        w(col)*((shift_pre(r,c,col) - ...
                                        data.T(row_index,col_index,col))^2);
                                end
                            end
                        end
                    end
                    
                    if (current_error < lowest_error),
                        lowest_error = current_error;
                        best_row = row_t;
                        best_col = col_t;
                    end
                    
                end
                
                % update the source_map with newly chosen pixel
                source_map(ypos,xpos) = (best_row - 1)*data.cols_t + best_col;
                
                % bestypos, bestxpos - pixel coords of best match in T
                bestypos = best_row;
                bestxpos = best_col;

                % paste that pixel into out and pre
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

% private function. used frequently in overlap_resynthesis_knn.m
function wrapped = add_and_wrap(value, accum, max_value)

wrapped = value + accum;
if (wrapped > max_value), 
    wrapped = wrapped - max_value; 
end
if (wrapped < 1), 
    wrapped = max_value + wrapped; 
end
