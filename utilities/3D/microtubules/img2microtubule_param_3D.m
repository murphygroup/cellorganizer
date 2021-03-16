function [mt_param, options] = img2microtubule_param_3D(im_mt, im_dna_seg, im_cell_seg, savedir, options)
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

% Copyright (C) 2015-2016 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

warning(['THE FUNCTION ' mfilename ' IS INCOMPLETE. Please read file comments for details.'])

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

if ~isfield(options.model, 'microtubule')
    options.model.microtubule = true;
end

if ~isfield(options.model.microtubule, 'searchparams')
    options.model.microtubule.searchparams = true;
end

options.model.microtubule = ml_initparam(options.model.microtubule, ...
    struct( ...
    'useCurrentResults', false ...
    ));

options = ml_initparam(options, struct(...
    'verbose', false ...
    ));

options.model.microtubule.searchparams = ml_initparam(options.model.microtubule.searchparams, ...
    struct( ...
    'n', [5,50:50:350,400,450], ...
    'mulen', [5 10 15 20 25 30 35 40 45], ...
    'colli_min_number', collinearityrange(5,0.97), ...
    'coeff_var', 1, ... %[0, 0.1, 0.2, 0.3], ... %this parameter does not get used
    'img_framework', 1, ...
    'replicate', 1 ... %this is a vector, we just use these as seeds for the microtubule growing
    ));


options.mt_dir = savedir;
if ~exist(options.mt_dir, 'dir')
    mkdir(options.mt_dir)
end

options.model.microtubule.featdir =  [options.mt_dir filesep 'img_feats'];
if ~exist(options.model.microtubule.featdir, 'dir')
    mkdir(options.model.microtubule.featdir);
end

options.model.microtubule.featdir_synth = [options.mt_dir filesep 'img_feats_synth'];
if ~exist(options.model.microtubule.featdir_synth, 'dir')
    mkdir(options.model.microtubule.featdir_synth)
end

%Loop through and calculate image features for all of images
imagefeat_path = [options.model.microtubule.featdir filesep 'feats' ];

[can_start, imagefeat_path, file_exists, tmpfile] = chunk_start(imagefeat_path, '.mat');
if can_start && ~file_exists %grj 11/3/14, chunk_start this later
    
    centrosome_coord = get_centrosome_coord(im_mt, im_dna_seg, im_cell_seg);
    
    disp(['Computing image features.'])
    [imfeats,idxes] = img2microtubule_feats(im_mt,im_cell_seg, im_dna_seg, centrosome_coord);
    
    save(imagefeat_path, 'imfeats', 'centrosome_coord', 'idxes')
    
    chunk_finish(tmpfile)
elseif file_exists
    load(imagefeat_path)
end

searchparams = struct2cell(options.model.microtubule.searchparams);
n_params = cellfun(@length, searchparams)';

% XYZres = [0.2 0.2 0.2]; %this is a default value from the microtubule code
XYZres.cell = options.model.resolution;
XYZres.objects = options.model.resolution;
XYZres.centrosome = options.model.resolution;


%this should be fixed eventually - grj 5/4/15
warning('Using default MT learning PSF')
warning('Please consider learning from data');
PSF = get_mt_psf(); %%%THIS PSF IS DATASET SPECIFIC, LEARN IT FROM TEH DATA IF IT DOES NOT EXIST


param_tmp = options;
param_tmp.skipimg = true;

%randomly shuffle the images to synthesize, so when this is done in
%parallel, the threads dont "pile up" on the same jobs

rng('default')
rng('shuffle')

%Shuffle the indexes so the results are unbiased if we choose to stop
%modeling
shuffle_inds = randperm(prod(n_params));
for ind = 1:prod(n_params)
    i = shuffle_inds(ind);
    
    param_coord = ind2sub_vec(n_params,i);
    
    
    img_num = sub2ind(n_params, param_coord(1), param_coord(2), param_coord(3), param_coord(4), param_coord(5), param_coord(6));
    
    imfeats_synth_files{ind} = [options.model.microtubule.featdir_synth filesep 'img' num2str(img_num)];
    
    [can_start, imfeats_synth_files{ind}, final_exists, tmpfile] = chunk_start(imfeats_synth_files{ind}, '.mat');
    
    if can_start & ~final_exists
        param_coord = ind2sub_vec(n_params, i);
        
        n =                 searchparams{1}(param_coord(1));
        mu_len =            searchparams{2}(param_coord(2));
        colli_min_number =  searchparams{3}(param_coord(3));
        coeff_var =         searchparams{4}(param_coord(4));
        img_framework =     searchparams{5}(param_coord(5));
        replicate =         searchparams{6}(param_coord(6)); %we use the replicate number as the random seed
        
        
        %             centrosome_coord = [];
        %             load(imagefeat_paths{img_framework}, 'centrosome_coord')
        
        disp(['Synthesizing images for ' num2str([param_coord]) '.']);
        
        valid_region = im_cell_seg & ~im_dna_seg;
        
        
        [im_mt_synth,imXYZ,mtXYZ,randlengths,resolution] = MT_synth(n, mu_len,colli_min_number, ~valid_region,centrosome_coord,XYZres,'bounce_control', replicate);
        %                     [im_mt,tcheck,imXYZ,mtXYZ,randlengths] = colli_generator_acute_jl2(n, mu_len, sigma_len, colli_min_number, ~valid_region, centrosome_coord, XYZres, 0.3, XYZres(1), replicate ,0,0,0,2, 0, 0,0);
        %                     [im_mt,tcheck,imXYZ,mtXYZ,randlengths] = colli_generator_acute_jl2_2(n,mu_len,sigma_len,colli_min_number,~valid_region,centrosome_coord,XYZres,0.3,XYZres(1),replicate);
        
        if isempty(im_mt_synth)
            continue;
        end
        im_mt_psf = imfilter(im_mt_synth, PSF);
        
        %Get features for that pattern
        imfeats_synth = img2microtubule_feats(im_mt_psf, im_cell_seg, im_dna_seg, centrosome_coord);
        
        %here we don't save the image to save storage
        save(imfeats_synth_files{ind}, 'imfeats_synth', 'n', 'mu_len', 'coeff_var', 'colli_min_number', 'img_framework', 'replicate', 'centrosome_coord', 'param_coord') %maybe also save the image
        
        chunk_finish(tmpfile)
    end
    %          end
end

%
% imfeats_real = [];
% c = 1;
%
% for i = 1:length(image_paths)
%     imfeats_file = [options.model.microtubule.featdir filesep 'img' num2str(i) '.mat'];
%     if exist(imfeats_file, 'file')
%         load(imfeats_file)
%         %load the file into the feature matrix
%         imfeats_real(c,:) = imfeats;
%         imfeats_imnum(c) = i;
%         c = c+1;
%     end
% end

imfeats_real = load(imagefeat_path, 'imfeats');
imfeats_real = imfeats_real.imfeats;

imfeats_real(isnan(imfeats_real)) = 0;

%Because we're dealing with a lot of image features, we first PCA the real
%image features, and project all the synthesized features into that PCA
%space. We use the PCA features that account for 95% of the variance. This
%is NOT in the original paper, but should work reasonably well with smaller
%datasets. GRJ 11/11/14

% [z, mu, sigma] = zscore(imfeats_real);
%
% [coeff, pcafeats, latent] = pca(z, 'Centered', false);
%
% ind = find((cumsum(latent)/sum(latent)) > 0.95);
% pcaind = ind(1);
%
% pcafeats = pcafeats(:, 1:pcaind);
%
% cov_sigma = cov(pcafeats(:,1:pcaind));


%%%this is not correctly done. Must match image to nearest synthesized
%%%image

c = 1;

imfeats_synth_all = [];
coordinate_inds = [];
dists = [];
%for each parameter coordinate

for i = 1:length(imfeats_synth_files)
    if exist(imfeats_synth_files{i}, 'file')
        synth_param(c) = load(imfeats_synth_files{i});
        c= c+1;
    end
end

imfeats_synth = vertcat(synth_param.imfeats_synth);

feats = zscore([imfeats_real;imfeats_synth]);

imfeats_real = feats(1,:);
imfeats_synth = feats(2:end,:);

[~,min_ind] = min(pdist2(imfeats_real, imfeats_synth));

mt_param.type = 'network';
mt_param.class = 'microtubule';
mt_param.n = synth_param(min_ind).n;
mt_param.mu_len = synth_param(min_ind).mu_len;
mt_param.colli_min_number = synth_param(min_ind).colli_min_number;
mt_param.version = '1.2';

end

function centrosome_coordinate = get_centrosome_coord(im_mt, im_dna_seg, im_cell_seg)
%this will have to be changed if the preprocessing folder structure changes
%obviously this is not the correct way to do it...
if ~exist('im_dna_seg', 'var')
    im_dna_seg = false(size(im_mt));
end

if ~exist('im_cell_seg', 'var')
    im_cell_seg = true(size(immt));
end

centrosome_coordinate = img2centrosome_coord(im_mt);

%if the centrosome coordinate is outside of the valid region,
%return the nearest point to the centrosome

valid_region = ~im_dna_seg & im_cell_seg;
if valid_region(centrosome_coordinate(1), centrosome_coordinate(2), centrosome_coordinate(3)) == 0
    [x,y,z] = ind2sub(size(valid_region), find(valid_region));
    [~, ind] = min(pdist2(centrosome_coordinate, [x,y,z]));
    centrosome_coordinate = [x(ind), y(ind), z(ind)];
end
end