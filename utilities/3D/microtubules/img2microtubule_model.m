function [model, param] = img2microtubule_model(image_paths, param)
% Implementation of microtubule learning model from 
%
%   J. Li, A. Shariff, M. Wiking, E. Lundberg, G.K. Rohde and R.F. Murphy 
%   (2012) Estimating microtubule distributions from 2D immunofluorescence 
%   microscopy images reveals differences among human cultured cell lines. 
%   PLoS ONE 7:e0050292.
%
% Adapted from a combination of software from
%
%   http://murphylab.web.cmu.edu/software/2012_PLoS_ONE_Microtubule_Models
%
% and updated models and functions currently in CellOrganizer 
%
%Gregory Johnson 11/9/14

warning(['THE FUNCTION ' mfilename 'IS INCOMPLETE. Please file comments for details'])

%%% Current things that are NOT implemented:
%       Centrosome model learning.
%       Not sure if any x,y,z scaling is implemented.
%       Image intensity adjustment to that learned from MTs.
%       Centrosome finding is currently the brightest spot in the cell, NOT
%         what is described in the paper.
%       Some parameters are hardcoded, should be impelmented in param
%         structure.
%       
%%% Known Bugs
%       Sometimes the image synthesis will fail. This is likely due to the
%         lack of 'space' in the cytoplasm in some segmented images.
%       Covariance of real image features sometimes results in an
%          ill-conditioned covariance matrix. Here we use 95% of Principal 
%          components to reduce the number of features.
%
%   As of 11/11/14 This code has not be tested on anything other than the
%   Tub images from the HeLa dataset in demo3Dimg2microtubule_model
%
%Gregory Johnson 11/11/14

warning off

if ~isfield(param.model, 'microtubule') 
    param.model.microtubule = true;
end

if ~isfield(param.model.microtubule, 'searchparams')
    param.model.microtubule.searchparams = true;
end

param.model.microtubule = ml_initparam(param.model.microtubule, ...
    struct( ...
    'useCurrentResults', false ...
    ));

param = ml_initparam(param, struct(...
    'verbose', false ...
    ));

param.model.microtubule.searchparams = ml_initparam(param.model.microtubule.searchparams, ...
    struct( ...
    'n', [5,50:50:350,400,450], ...
    'mulen', [5 10 15 20 25 30 35 40 45], ...
    'colli_min_number', collinearityrange(5,0.97), ...
    'coeff_var', 1, ... %[0, 0.1, 0.2, 0.3], ... %this parameter does not get used
    'img_framework', 1:length(image_paths), ...
    'replicate', 1 ... %this is a vector, we just use these as seeds for the microtubule growing
    ));


param.mt_dir = [param.tempparent filesep 'microtubule'];
if ~exist(param.mt_dir, 'dir')
    mkdir(param.mt_dir)
end

param.model.microtubule.featdir =  [param.mt_dir filesep 'img_feats'];
if ~exist(param.model.microtubule.featdir, 'dir')
    mkdir(param.model.microtubule.featdir);
end

param.model.microtubule.featdir_synth = [param.mt_dir filesep 'img_feats_synth'];
if ~exist(param.model.microtubule.featdir_synth, 'dir')
    mkdir(param.model.microtubule.featdir_synth)
end


imagefeat_paths = cell(1, length(image_paths));

%Loop through and calculate image features for all of images
for i = 1:length(image_paths) %this can be a chunk_start loop and/or function, This can be in teh preprocessing loop too
    imagefeat_paths{i} = [param.model.microtubule.featdir filesep 'img' num2str(i)];

    [can_start, imagefeat_paths{i}, file_exists, tmpfile] = chunk_start(imagefeat_paths{i}, '.mat');
    if can_start & ~file_exists %grj 11/3/14, chunk_start this later
        
        [improt, cellmask, dnamask, centrosome_coord] = get_seg_prot(image_paths{i},i, param);
        
        if isempty(improt)
            continue;
        end
        
        disp(['Computing image features for image ' num2str(i) ' of ' num2str(length(image_paths))])
        [imfeats,idxes] = img2microtubule_feats(improt,cellmask, dnamask, centrosome_coord);
        
        save(imagefeat_paths{i}, 'imfeats', 'centrosome_coord', 'idxes')
        
        chunk_finish(tmpfile)
    end 
end

searchparams = struct2cell(param.model.microtubule.searchparams);
n_params = cellfun(@length, searchparams)';

% XYZres = [0.2 0.2 0.2]; %this is a default value from the microtubule code
XYZres.cell = param.model.resolution;
XYZres.objects = param.model.resolution;
XYZres.centrosome = param.model.resolution;


%this should be fixed eventually - grj 5/4/15
warning('Using default MT learning PSF')
warning('Please consider learning from data');
PSF = get_mt_psf(); %%%THIS PSF IS DATASET SPECIFIC, LEARN IT FROM TEH DATA IF IT DOES NOT EXIST


param_tmp = param;
param_tmp.skipimg = true;

%randomly shuffle the images to synthesize, so when this is done in
%parallel, the threads dont "pile up" on the same jobs

if ~param.model.microtubule.useCurrentResults
    
    rng('default')
    rng('shuffle')
    
    
    %Shuffle the indexes so the results are unbiased if we choose to stop
    %modeling
    shuffle_inds = randperm(prod(n_params));
    for ind = 1:prod(n_params)
        i = shuffle_inds(ind);

        param_coord = ind2sub_vec(n_params,i);


        img_num = sub2ind(n_params, param_coord(1), param_coord(2), param_coord(3), param_coord(4), param_coord(5), param_coord(6));

        imfeats_synth_file = [param.model.microtubule.featdir_synth filesep 'img' num2str(img_num)];

        if ~exist([imfeats_synth_file '.mat'], 'file')
            [can_start, imfeats_synth_file, final_exists, tmpfile] = chunk_start(imfeats_synth_file, '.mat');

            if can_start & ~final_exists
                param_coord = ind2sub_vec(n_params, i);

                n =                 searchparams{1}(param_coord(1));
                mu_len =            searchparams{2}(param_coord(2));
                colli_min_number =  searchparams{3}(param_coord(3));
                coeff_var =         searchparams{4}(param_coord(4));
                img_framework =     searchparams{5}(param_coord(5));
                replicate =         searchparams{6}(param_coord(6)); %we use the replicate number as the random seed

                [~, cellmask, dnamask] = get_seg_prot(image_paths{img_framework},img_framework, param_tmp);
                if isempty(cellmask)
                    continue
                end

                centrosome_coord = [];
                load(imagefeat_paths{img_framework}, 'centrosome_coord')

                disp(['Synthesizing images for ' num2str([param_coord]) '.']);

                valid_region = cellmask & ~dnamask;
                
                
                [im_mt,imXYZ,mtXYZ,randlengths,resolution] = MT_synth(n, mu_len,colli_min_number, ~valid_region,centrosome_coord,XYZres,'bounce_control', replicate);
%                     [im_mt,tcheck,imXYZ,mtXYZ,randlengths] = colli_generator_acute_jl2(n, mu_len, sigma_len, colli_min_number, ~valid_region, centrosome_coord, XYZres, 0.3, XYZres(1), replicate ,0,0,0,2, 0, 0,0);
%                     [im_mt,tcheck,imXYZ,mtXYZ,randlengths] = colli_generator_acute_jl2_2(n,mu_len,sigma_len,colli_min_number,~valid_region,centrosome_coord,XYZres,0.3,XYZres(1),replicate);

                if isempty(im_mt)
                    continue;
                end
                im_mt_psf = imfilter(im_mt, PSF);

                %Get features for that pattern
                imfeats_synth = img2microtubule_feats(im_mt_psf, cellmask, dnamask, centrosome_coord);

                %here we don't save the image to save storage
                save(imfeats_synth_file, 'imfeats_synth', 'n', 'mu_len', 'coeff_var', 'colli_min_number', 'img_framework', 'replicate', 'centrosome_coord', 'param_coord') %maybe also save the image

                chunk_finish(tmpfile)
            end
        end
%          end
    end
end


imfeats_real = [];
c = 1;

for i = 1:length(image_paths)
    imfeats_file = [param.model.microtubule.featdir filesep 'img' num2str(i) '.mat'];
    if exist(imfeats_file, 'file')
        load(imfeats_file)
        %load the file into the feature matrix
        imfeats_real(c,:) = imfeats;
        imfeats_imnum(c) = i;
        c = c+1;
    end
end

imfeats_real(isnan(imfeats_real)) = 0;


%Because we're dealing with a lot of image features, we first PCA the real
%image features, and project all the synthesized features into that PCA
%space. We use the PCA features that account for 95% of the variance. This
%is NOT in the original paper, but should work reasonably well with smaller
%datasets. GRJ 11/11/14

[z, mu, sigma] = zscore(imfeats_real);

[coeff, pcafeats, latent] = pca(z, 'Centered', false);

ind = find((cumsum(latent)/sum(latent)) > 0.95);
pcaind = ind(1);

pcafeats = pcafeats(:, 1:pcaind);

cov_sigma = cov(pcafeats(:,1:pcaind));


%%%this is not correctly done. Must match image to nearest synthesized
%%%image

%for every cell
q = 1;
for i = 1:length(imfeats_imnum)
    if param.verbose
        disp(['Finding nearest match for cell ' num2str(i) filesep num2str(length(imfeats_imnum))]);
    end
    
    param_cellnum = imfeats_imnum(i);
    
    c = 1;
    
    imfeats_synth_all = [];
    coordinate_inds = [];
    dists = [];
    %for each parameter coordinate
    for j = 1:prod(n_params(1:4))
        for k = 1:n_params(6)

        param_coord = ind2sub_vec(n_params(1:4),j);
        
        param_rep = k;
        
        coordinate_inds(c) = sub2ind(n_params, param_coord(1), param_coord(2), param_coord(3), param_coord(4), param_cellnum, param_rep);
        
        imfeats_synth_file = [param.model.microtubule.featdir_synth filesep 'img' num2str(coordinate_inds(c)) '.mat'];
        
        
            if exist(imfeats_synth_file, 'file')
                load(imfeats_synth_file)

                imfeats_synth_pca = ((imfeats_synth - mu)./sigma)*coeff;
                imfeats_synth_pca = imfeats_synth_pca(1:pcaind);
                
                imfeats_synth_all(c,:) = imfeats_synth;
                dists(c) = sqrt((pcafeats(i,:) - imfeats_synth_pca)*inv(cov_sigma)*(pcafeats(i,:) - imfeats_synth_pca)');
                
                c = c+1;
            end
        end
    end
    
    if ~isempty(dists)
        [~, min_ind] = min(dists);

        coordinate_vec= ind2sub_vec(n_params, coordinate_inds(min_ind));
        coord_final(q,:) = [searchparams{1}(coordinate_vec(1)), searchparams{2}(coordinate_vec(2)), searchparams{3}(coordinate_vec(3))];
        cell_ind(q) = param_cellnum;
        q = q+1;
    end
end
        
param.model.microtubule.cell_param = coord_final;
param.model.microtubule.cell_ind = cell_ind;

model.type = 'network';
model.class = 'microtubule';
model.resolution = param.model.resolution;
model.version = '1.1';
model.parameters.mean = mean(coord_final);
model.parameters.cov = cov(coord_final);

end

function [improt_seg, imcell, imdna, centrosome_coordinate] = get_seg_prot(impath, imnum, param)
%this will have to be changed if the preprocessing folder structure changes
%obviously this is not the correct way to do it...
    param = ml_initparam(param, struct('skipimg', false));
    
    
    segfile = [param.preprocessingFolder filesep 'cell' num2str(imnum) '.mat'];
    
    if ~exist(segfile, 'file')
        improt_seg = [];
        imcell = [];
        imdna = [];
        centrosome_coordinate = [];
        return
    end
    
    load(segfile)
    
    if all(segcell(:) == 0) || all(segdna(:) == 0)
        improt_seg = [];
        imcell = [];
        imdna = [];
        centrosome_coordinate = [];
        return
    end
    
    imcell = zeros(original_segcellsize);
    imcell(:,:,bot_slice:top_slice) = segcell;
    
    imdna = zeros(original_segcellsize);
    imdna(:,:,bot_slice:top_slice) = segdna;
    
    if ~param.skipimg
        improt_res = ml_downsize(ml_readimage(impath), param.downsampling);
        
        improt_seg = zeros(original_segcellsize);
        improt_seg(:,:,bot_slice:top_slice) = segcell.*improt_res(:,:,bot_slice:top_slice);
    
        centrosome_coordinate = img2centrosome_coord(improt_seg);
        
        %if the centrosome coordinate is outside of the valid region,
        %return the nearest point to the centrosome
        valid_region = ~segdna & segcell;
        if valid_region(centrosome_coordinate(1), centrosome_coordinate(2), centrosome_coordinate(3)) == 0
            [x,y,z] = ind2sub(size(valid_region), find(valid_region));
            [~, ind] = min(pdist2(centrosome_coordinate, [x,y,z]));
            centrosome_coordinate = [x(ind), y(ind), z(ind)];
        end
        
    else
        improt_seg = [];
    end
end




