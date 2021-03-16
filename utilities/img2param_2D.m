function param = img2param_2D( imdna_path,imcell_path,...
    improt_path,immask_path, savedir, options )

% Gregory Johnson
%
% Copyright (C) 2016 Murphy Lab
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

% January 13, 2018 I. Cao-Berg Cleaned up method
% Feb 17, 2018 X. Ruan Add 2d PCA model

options = ml_initparam(options, struct('downsampling', [1,1], ...
    'display', false, ...
    'debug', false, ...
    'verbose', false, ...
    'train', []));

options.train = ml_initparam(options.train, struct('flag', 'all'));

options.model.resolution_orig = options.model.resolution;
options.model.resolution = options.model.resolution.*options.downsampling;

options = set_default_options(options);

disp('Loading images');
imdna = loadImage(imdna_path, options.downsampling);
imcell = loadImage(imcell_path, options.downsampling);
improt = loadImage(improt_path, options.downsampling);
immask = loadImage(immask_path, options.downsampling);

param.imsize = size(imdna);

disp('Segmenting image');
savefile = [savedir filesep 'seg.mat'];
if ~exist(savefile, 'file')
    contrast_stretch_image = @(x)(imadjust(x,stretchlim(x),[]));
    
    is_it_contrast_stretchable = @(x)(strcmpi(class(x),'uint8') && ...
        numel(unique(x)) == 2 );
    
    if is_it_contrast_stretchable(imdna)
        disp('Contrast stretching nuclear membrane image');
        imdna = contrast_stretch_image( imdna );
    end
    
    if is_it_contrast_stretchable(imcell)
        disp('Contrast stretching cell membrane image');
        imcell = contrast_stretch_image( imcell );
    end
    
    [seg_dna, seg_cell] = seg_cell_and_nucleus_2D( ...
        imdna, ...
        imcell, ...
        immask, ...
        options);
    
    disp('Saving segmentations to disk');
    is_correct = @(x)(strcmpi(class(x),'double') && ...
        numel(unique(x)) == 2 && max(max(x)) ~= 1 );
    
    if ~is_correct( seg_dna )
        disp('Segmentation should be saved as double. Casting to double.');
        seg_dna = double(seg_dna);
        seg_dna(find(seg_dna~=0)) = 1;
        seg.nuc = seg_dna;
    else
        seg.nuc = seg_dna;
    end
    
    if ~is_correct( seg_cell )
        disp('Segmentation should be saved as double. Casting to double.');
        seg_cell = double(seg_cell);
        seg_cell(find(seg_cell~=0)) = 1;
        seg.cell = seg_cell;
    else
        seg.cell = seg_cell;
    end
    
    save(savefile, 'seg')
else
    load(savefile)
end
param.seg = seg;

%Compute Nuclear image feats
% xruan 01/05/2016 add the case for diffeomorphic
savefile = [savedir filesep 'nuc.mat'];
if ~isempty(imdna)
    
    options.nucleus = ml_initparam(options.nucleus, struct('type', 'medial axis'));
    
    if exist(savefile, 'file')
        load(savefile)
    else
        switch options.nucleus.type
            case 'medial_axis'
                nuc = img2medial_axis(param.seg.nuc);
            case 'pca'
                nuc = img2landmarks(param.seg.nuc);                
            case 'diffeomorphic'
                nuc.type = 'diffeomorphic';
                nuc.seg = param.seg.nuc;
            otherwise
                error(['Unrecognized nuclear shape model type ' options.protein.type])
        end
        
        save(savefile, 'nuc');
    end
else
    load(savefile)
end
% xruan 03/08/2016 change nucleus to nucleus to unify the names with 3D
% case
param.nuc = nuc;

%Now get the Cell level parameters
% xruan 01/05/2016 add the case for diffeomorphic
savefile = [savedir filesep 'cell.mat'];
if ~isempty(imcell)
    
    options.cell = ml_initparam(options.cell, struct('type', 'ratio'));
    
    if ~exist(savefile, 'file')
        switch options.cell.type
            case 'ratio'
                options.nuclear_image_filename = imdna_path;
                options.cell_image_filename = imcell_path;
                cellfit = ml_parsecell2D({}, ...
                    ml_findmainobj2d_bw(param.seg.cell), ...
                    ml_findmainobj2d_bw(param.seg.nuc), ... 
                    1,param.imsize,...
                    {'da','nucarea','nuccenter','nucmangle','nuchitpts',...
                    'nuccontour','nucellhitpts','nucdist','nucelldist','nucecc',...
                    'cellarea','cellcenter','cellmangle','cellcontour',...
                    'cellhitpts','celldist','cellecc'}, options );
            case 'pca'
                cellfit = img2landmarks(param.seg.cell);                                
            case 'diffeomorphic'
                cellfit.type = 'diffeomorphic';
                cellfit.seg = param.seg.cell;
            otherwise
                error(['Unrecognized nuclear shape model type ' options.protein.type])
        end
        save(savefile, 'cellfit')
    else
        load(savefile)
    end
    param.cell = cellfit;
end

savefile = [savedir filesep 'preprocess.mat'];
tic
if ~isempty(improt) && ~exist(savefile, 'file')
    if options.debug && options.display
        imshow( improt );
        title( imdna )
        %interactive_display(options);
    end
    
    options = ml_initparam(options, struct(...
        'threshmeth_prot', 'nih'));
    procimage = ...
        ml_preprocess(double(improt),immask,'ml','yesbgsub',options.threshmeth_prot);
    
    save(savefile, 'procimage')
elseif exist(savefile, 'file')
    load(savefile)
else
    procimage = [];
end
param.prot.procimage = procimage;

%Now compute the protein models.
if ~isempty(improt)
    savefile = [savedir filesep 'prot.mat'];
    
    options.protein = ml_initparam(options.protein, struct('type', 'vesicle'));
    
    if exist(savefile, 'file')
        load(savefile)
    else
        switch options.protein.type
            case 'gmm'
                if ~isfield(options, 'ml_objgaussmix')
                    options.ml_objgaussmix = [];
                end
                if options.display
                    options.ml_objgaussmix.isshow = true;
                end
                prot = img2protfit_2D(procimage, options.ml_objgaussmix);
            case 'network'
                %%% do this for the microtubule model %%%
            otherwise
                error(['Unrecognized protein model type ' options.protein.type])
        end
        save(savefile, 'prot')
    end
    param.prot = prot;
end


param.options = options;
end

function options = set_default_options(options)
options = ml_initparam(options, struct('cell', [], ...
    'nucleus', [], ...
    'protein', []));

if ~isfield(options.nucleus, 'type') || isempty(options.nucleus.type)
    options.nucleus.type = 'medial axis';
end

if ~isfield(options.cell, 'type') || isempty(options.cell.type)
    options.cell.type = 'ratio';
end

if ~isfield(options.protein, 'type') || isempty(options.protein.type)
    options.protein.type = 'vesicle';
end
end