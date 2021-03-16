function [struct_euc, struct_orig, param, positions] = compute_embedding_dimensionality( distance_matrix, param )
%Finds the "best" best embedding dimension of a distance matrix.
%
%Computes an embedding at multiple dimensionalities, and measures the 1-r^2
%compared to the original distance matrix, and compared to a euclidian
%version of that matrix, and returns the dimension that is one above the
%the dimension with the largest positive 2nd derivative compared to the
%residuals with greater and lower dimensionality.
%
%The 2nd derivative currently does not account for skipping dimensions

%gregory r johnson, 2/28/15

if ~exist('param', 'var') | isempty(param);
    param = struct();
end

param = ml_initparam(param, struct('method',  'regress', ...
                                    'weight_factor', 0, ...
                                    'embedding_dims', [1:15], ...
                                    'paralleldir', [], ...
                                    'calceuc', true, ...
                                    'skip_missing', false));

if ~isempty(param.paralleldir)
    mkdir(param.paralleldir)
end

%make sure the dimensions are in ascending order
param.embedding_dims = sort(param.embedding_dims);

if param.calceuc
    param.desired_shape_space_dimensionality = inf;
    pos_euc = embed_partial_distance_matrix(distance_matrix, param);
    dists_euc = squareform(pdist(pos_euc));
end

resid_euc = nan(1, length(param.embedding_dims));
resid_orig = nan(1, length(param.embedding_dims));

for i = 1:length(param.embedding_dims)
    
    param.desired_shape_space_dimensionality = param.embedding_dims(i);
    
    %if there is a parallelsave directory
    %do this very convoluted thing
    if ~isempty(param.paralleldir)
        clear positions
        clear dists
        
        savefile = [param.paralleldir filesep 'embed_' num2str(i)];

        [can_start, savefile, final_exists, temp_name] = chunk_start( savefile);
        if can_start & ~final_exists & ~param.skip_missing
            [pos_temp, ~, ~] = embed_partial_distance_matrix(distance_matrix, param);    
            
            save(savefile, 'pos_temp')

            chunk_finish(temp_name)

            positions{i} = pos_temp;
            dists{i} = squareform(pdist(pos_temp));
        elseif final_exists
            load(savefile)
            positions{i} = pos_temp;
            dists{i} = squareform(pdist(pos_temp));
        end

    %otherwise just do it all in one thread
    else    
        [positions{i}, ~, ~] = embed_partial_distance_matrix(distance_matrix, param);    
        dists{i} = squareform(pdist(positions{i}));
    end
    
    if exist('positions', 'var')
    
        nans = isnan(distance_matrix);

        if param.calceuc
           r = corrcoef(dists{i}(~nans), dists_euc(~nans));
           r = r(1,2);
           resid_euc(i) = 1 - r^2;
        end
        
        r = corrcoef(dists{i}(~nans), distance_matrix(~nans));
        r = r(1,2);

        resid_orig(i) = 1 - r^2;
    end
end


%if any of the computations are missing
if any(isnan(resid_orig)) && ~param.skip_missing 
    struct_euc = [];
    struct_orig = [];
    return
else

    if param.calceuc
        del2_euc =  gradient(gradient(resid_euc));
        %exclude the first coordinate
        [~, ind] = max(del2_euc(2:end));
        ind_euc = ind+2;
        dim_euc = param.embedding_dims(ind_euc);

        struct_euc.dim = dim_euc;
        struct_euc.resid = resid_euc;
        struct_euc.del2 = del2_euc;
        
        if ~isempty(param.paralleldir)
            load([param.paralleldir filesep 'embed_' num2str(ind_euc) '.mat'])
            struct_euc.pos = pos_temp;
        else
            struct_euc.pos = positions{ind_euc};
        end
        
        struct_euc.pos = positions{ind_euc};
        struct_euc.pos_best = pos_euc;
    else
        struct_euc = [];
    end


    del2_orig =  gradient(gradient(resid_orig(2:end)));
    %exclude the first coordinate
    [~, ind] = max(del2_orig(2:end));
    ind_orig = ind+2;
    dim_orig = param.embedding_dims(ind_orig);

    struct_orig.dim = dim_orig;
    struct_orig.resid = resid_orig;
    struct_orig.del2 = del2_orig;
    if ~isempty(param.paralleldir)
        load([param.paralleldir filesep 'embed_' num2str(ind_orig) '.mat'])
        struct_orig.pos = pos_temp;
    else
        struct_orig.pos = positions{ind_orig};
    end
    struct_orig.pos_best = distance_matrix;
end


end
