function [object,resolution] = model2tcell(nucleus, cellMembrane, model,param)
% The main function for t cell protein pattern synthesis. 
% Input: nucleus: nucleus shape; cellMembrane: cell shape; model: protein model for t cell; param: input options.
% Output: object: t cell protein pattern,resolution: resolution of the model. 
% 
% It contains the pipeline for the synthesis. 
% tcell_voxel_intensity_sampling.m is used for sampling protein pattern from the t cell model. 
% tcell_cell_shape_matching.m is used for mapping the protein image in the template to the input cell shape. 
%
%
% Author: Xiongtao Ruan 
% Date: Sept. 17, 2016

% sampling the intensities for the voxels
[I_raw, I_seg, param] = tcell_voxel_intensity_sampling(model, param);

% matching the shape to the input shape
[I_match, param] = tcell_cell_shape_matching(I_seg, I_raw, cellMembrane, model, param);

resolution = param.resolution.objects; 
object = I_match;


end