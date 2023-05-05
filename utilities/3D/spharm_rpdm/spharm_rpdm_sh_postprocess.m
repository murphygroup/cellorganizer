function [cell_post, nuc_post] = spharm_rpdm_sh_postprocess(cellfit, nuc, savedir, options) 
% postprocess function for both cell and nuclear parameterizations. 
% 
% 09/19/2018 xruan: add code for filtering bad parameterization, and
% postprocess for one component.
% 02/24/2019 xruan: pass alignment option for postprocess
% 02/24/2019 xruan: update postprocess setting in terms of train flag
% 03/23/2023 R.F.Murphy pass final hausdorff distance and jaccard index to postprocessing
% 04/12/2023 R.F.Murphy fix error in calling sequence introduced by previous change
% 5/4/2023 R.F. Murphy also fix calling sequence for nuclear models

if strcmp(options.train.flag, 'cell') || strcmp(options.train.flag, 'framework') || strcmp(options.train.flag, 'all')
    cell_postprocess_savefile = [savedir filesep 'cell_postprocess.mat'];
    if ~exist(cell_postprocess_savefile, 'file')
        % postprocess_options = struct();
        postprocess_options = options.spharm_rpdm;
        [cell_post] = spherical_parameterization_postprocess(cellfit.vertices, cellfit.faces, cellfit.sph_verts, cellfit.fvec, cellfit.final_hd, cellfit.jaccard_index, postprocess_options);
        save(cell_postprocess_savefile, 'cell_post');
    else
        load(cell_postprocess_savefile, 'cell_post');
    end
    if  strcmp(options.train.flag, 'cell')
        nuc_post = [];
    end
end

if strcmp(options.train.flag, 'nuclear') || strcmp(options.train.flag, 'framework') || strcmp(options.train.flag, 'all')
    nuc_postprocess_savefile = [savedir filesep 'nuc_postprocess.mat'];        
    if ~exist(nuc_postprocess_savefile, 'file')
        if strcmp(options.train.flag, 'nuclear')
            postprocess_options = options.spharm_rpdm;
            [nuc_post] = spherical_parameterization_postprocess(nuc.vertices, nuc.faces, nuc.sph_verts, nuc.fvec, nucfit.final_hd, nucfit.jaccard_index, postprocess_options);
            save(nuc_postprocess_savefile, 'nuc_post');
        elseif strcmp(options.train.flag, 'framework') || strcmp(options.train.flag, 'all')
            if isempty(cell_post)
                nuc_post = [];
            else
                postprocess_options.use_given_rotation_matrix = true;
                postprocess_options.use_given_center = true;
                postprocess_options.rotation_matrix = cell_post.R;
                postprocess_options.rotation_center = cell_post.rotation_center;
                [nuc_post] = spherical_parameterization_postprocess(nuc.vertices, nuc.faces, nuc.sph_verts, nuc.fvec, nucfit.final_hd, nucfit.jaccard_index, postprocess_options);
            end
        end
        save(nuc_postprocess_savefile, 'nuc_post');
    else
        load(nuc_postprocess_savefile, 'nuc_post');
    end
    if  strcmp(options.train.flag, 'nuclear')
        cell_post = [];
    end    
end     
        
end
