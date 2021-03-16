function param = img2param_3D( imdna_path,imcell_path,...
    improt_path,immask_path, savedir, options )
% IMG2PARAM_3D Parameterizes an image of a cell for use in the CellOrganizer suite.
%
% Inputs
% ------
% imdna_path    cell array of image files for dna parameter extraction
% imcell_path   cell array of image files for cell parameter extraction
%               (optional if options.train.flag = 'nuclear')
% improt_path   cell array of image files for protein parameter extraction
%               (optional if options.train.flag = 'nuclear' or 'framework')
% immask_path   cell array of image files for masking a single cell
%               from a given image (these are always optional and can
%               be passed in as an array of [])
% options       a structure holding possible parameter options. This
%               struct should also contain all temporary results
%               folders.(see set_temp_result_folders.m)

% Author: Devin Sullivan 6/13/13
%
% Copyright (C) 2007-2019 Murphy Lab
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
% For additional information visit http://murphylab.web.cmu.edu/ or
% send email to murphy@cmu.edu

% Gregory Johnson 6/17/13
%   Modified for single image
%   Shortened some variable names
%   Formatted for readability
%
% Xiongtao Ruan 01/06/2016
%   change seg_cell_and_nucleus_2D to seg_cell_and_nucleus_3D.
%   add the diffeomorphic case of parameter setting.
%
% Ivan E. Cao-Berg 04/15/2016
%   Fixed bug where image being downsampled twice.
%
% R.F. Murphy 3/11/2017
%   Allow nuc model building even if imdna is empty
%   (which allows nuclear hole finding to work)
%
% Xiongtao Ruan 09/14/2018
%   add SPHARM_RPDM model parameterization
%   fix bug for volume calculation
%
% 12/8/2020 R.F.Murphy use options.spharm_rpdm for both nuc and cell
% 2/7/2021 R.F. Murphy don't skip spharm postprocessing if_skip_nuclear_cell
% 2/8/2021 R.F. Murphy fix NaN handling after parameterization; check hd_thresh

options = ml_initparam(options, struct('downsampling', [1,1,1], ...
    'display', false, ...
    'train', [], ...
    'segminnucfraction', 0.17));

options.train = ml_initparam(options.train, struct('flag', 'all'));

if options.verbose
    disp(['Loading images with downsampling set to: [' ...
        num2str(options.downsampling) ']']);
end
imdna = loadImage(imdna_path, options.downsampling);
imcell = loadImage(imcell_path, options.downsampling);
improt = loadImage(improt_path, options.downsampling);
immask = loadImage(immask_path, options.downsampling);


% disp('testing sparm_obj_percell_3D');
% param = sparm_obj_percell_3D(imdna_path,...
%                     imcell_path,...
%                     improt_path,...
%                     immask_path, ...
%                     options);
% 
% 
% disp('finished testing');return;

savefile = [savedir filesep 'seg.mat'];
if ~exist(savefile, 'file') && (~isfield(options,'if_skip_cell_nuclear_model') || ~options.if_skip_cell_nuclear_model)
    % xruan 01/06/2016
    % change seg_cell_and_nucleus_2D to seg_cell_and_nucleus_3D
    %    [seg_dna, seg_cell] = seg_cell_and_nucleus_2D( ...
    tic;
    %options.display=1;
    [seg_dna, seg_cell] = seg_cell_and_nucleus_3D( ...
        imdna, ...
        imcell, ...
        improt, ...
        immask, ...
        [1 1 1], ...   %% xruan 01/06/2015 this is just set to be compatible with other cases.
        options);
    options.running_time.segmentation = toc;
    
    seg.nuc = seg_dna;
    seg.cell = seg_cell;
    
    % xruan 09/14/2018
    disp('Computing nuclear and cell volume');
    cellvolume = double(sum(seg_cell(:)));
    nuclearvolume = double(sum(seg_dna(:)));
    ncratio = nuclearvolume/cellvolume;
    disp(['Nuclear/cell volume ratio is ' num2str(ncratio)]);
    if ncratio < options.segminnucfraction
        warning('Nuclear volume fraction too low.');
        param = [];
        return
    end
    
    save(savefile, 'seg')
elseif (isfield(options,'if_skip_cell_nuclear_model') && options.if_skip_cell_nuclear_model)
        seg.cell=[];
        seg.nuc=[];


    if exist('imdna')
        img.nuc = imdna;
    end
    
    if exist('imcell')
        img.cell = imcell;
    end
    
    if exist('improt')
        img.prot = improt;
    end
    
    save(savefile, 'seg','img')

else
    load(savefile)
end
param.seg = seg;

%Compute nuclear image feats
savefile = [savedir filesep 'nuc.mat'];
%if ~isempty(imdna) && ~exist(savefile, 'file')
if ~exist(savefile, 'file') && (~isfield(options,'if_skip_cell_nuclear_model') || ~options.if_skip_cell_nuclear_model)
    tic
    switch options.nucleus.type
        case 'cylindrical_surface'
            [spfeat,surfmap] = tp_nucimgfeat(param.seg.nuc, ...
                'cylsurf', options );
            
            nuc.spfeat = spfeat;
            nuc.surfmap = surfmap;
            % xruan 01/06/2016
        case 'diffeomorphic'
            nuc.type = 'diffeomorphic';
            nuc.seg = param.seg.nuc;
            
            % xruan 09/14/2018
        case 'spharm_rpdm'
            [nuc] = spharm_rpdm_image_parameterization(param.seg.nuc, options.spharm_rpdm);
            nuc.type = 'spharm_rpdm';
        otherwise
            warning(['Unsupported nuclear model type ' options.nucleus.type '. Returning empty model.'])
            nuc = [];
    end
    
    save(savefile, 'nuc')
    options.running_time.nuclear_image_parameterization = toc;
elseif (isfield(options,'if_skip_cell_nuclear_model') && options.if_skip_cell_nuclear_model)
    nuc = [];
else
    load(savefile)
end
param.nuc = nuc;

% xruan 03/08/2016 change cell to cellfit to unify name and resolve
% conflition with matlab built-in function
%Now get the Cell level parameters
savefile = [savedir filesep 'cell.mat'];
if ~isempty(imcell) && (~isfield(options,'if_skip_cell_nuclear_model') || ~options.if_skip_cell_nuclear_model)
    if ~exist(savefile, 'file')
        tic
        switch options.cell.type
            case 'ratio'
                
                cellfit = cellfit_percell( ...
                    param.seg.nuc,  ...
                    param.seg.cell, ...
                    options);
                % xruan 01/06/2016
            case 'diffeomorphic'
                cellfit.type = 'diffeomorphic';
                cellfit.seg = param.seg.cell;
                
                % xruan 09/14/2018
            case 'spharm_rpdm'
                if strcmp(options.train.flag, 'nuclear')
                    cellfit = nuc;
                else
                    [cellfit] = spharm_rpdm_image_parameterization(param.seg.cell, options.spharm_rpdm);
                    if cellfit.final_hd > options.hd_thresh
                        error('error of parameterization above hd_thresh')
                    end
%                    if (sum(isnan(cellfit.cost_mat))>1)
%                        disp('warning!cost_mat')
%                        assert(sum(isnan(cellfit.cost_mat))<1)
%                    end 
                    if (sum(isnan(cellfit.vertices))>0)
                        error('bad vertices')
                    end
                    if (sum(isnan(cellfit.sph_verts))>0)
                        error('bad sph_verts')
                    end
                    cellfit.type = 'spharm_rpdm';
                end
            otherwise
                warning(['Unsupported cell model type ' options.nucleus.type '. Returning empty model.'])
                cellfit = [];
        end
        
        save(savefile, 'cellfit')
        options.running_time.cell_image_parameterization = toc;
    else
        load(savefile)
    end
    param.cell = cellfit;
end

% 09/17/2018 add postprocess for spharm rpdm model.
% 02/24/2019 update the type setting so that it works for only cell/nuclear shape as spharm_rpdm
%if (strcmp(options.cell.type, 'spharm_rpdm') || strcmp(options.nucleus.type, 'spharm_rpdm')) && options.spharm_rpdm.postprocess && isfield(options,'if_skip_cell_nuclear_model') && ~options.if_skip_cell_nuclear_model
% 02/07/2021 fix the logic for if_skip
if (strcmp(options.cell.type, 'spharm_rpdm') || strcmp(options.nucleus.type, 'spharm_rpdm')) && options.spharm_rpdm.postprocess && (~isfield(options,'if_skip_cell_nuclear_model') || ~options.if_skip_cell_nuclear_model)
    [cell_post, nuc_post] = spharm_rpdm_sh_postprocess(cellfit, nuc, savedir, options);
    % set the parameter use postprocess parameters
    param.nuc = nuc_post;
    param.cell = cell_post;
end

%Now compute the protein models.
if ~isempty(improt) && (strcmpi( options.train.flag, 'all' )||strcmpi( options.train.flag, 'protein' ))
    savefile = [savedir filesep 'prot.mat'];
    if exist(savefile, 'file')
        load(savefile)
    else
        tic
        switch options.protein.type
            case 'spharm_obj' % build shape model and ppm (Shen Jin, Sep 2019)
                disp('Preprocessing image for SPHARM object model')
                prot = spharm_obj_percell_3D(imdna_path,...
                    imcell_path,...
                    improt_path,...
                    immask_path, ...
                    options,.....
                    param.seg);
            case 'gmm'
                prot = protfit_percell_3D(param.seg.nuc, ...
                    param.seg.cell, ...
                    improt, ...
                    param.seg.cell, ...
                    options);
            case 'microtubule_growth'
                savedir_mt = [savedir filesep 'mt_feats'];
                
                prot = img2microtubule_param_3D(improt, ...
                    param.seg.nuc, ...
                    param.seg.cell, ...
                    savedir_mt, ...
                    options);
            otherwise
                error(['Unrecognized protein model type ' options.protein.type])
        end
        save(savefile, 'prot')
        options.running_time.protein_image_parameterization = toc;
    end
    param.prot = prot;
end


savefile = [savedir filesep 'compartment.mat'];
%Lastly use any preprocessed results to create compartment
%models
if ~exist(savefile, 'file')
    compartment = compartment_percell(...
        param.seg.nuc, ...
        param.seg.cell,...
        improt, ...
        options);
    
    save(savefile, 'compartment')
else
    load(savefile)
end

param.compartment = compartment;
param.options = options;
end
