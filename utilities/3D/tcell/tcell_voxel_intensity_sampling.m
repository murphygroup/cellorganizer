function [I_raw, I_seg, param] = tcell_voxel_intensity_sampling(proteinModel, param)
% Input: proteinModel: the model for t cell; param: options for the pipeline.
% 
% The idea is to sampling the intensities for the voxels from the model. And we provide two methods:
% Gaussian: using the model mean and standard deviation to sample the voxels. 
% Empirical: using the empirical distribution of each voxel to sample the intensity. 
%
% Author: Xiongtao Ruan


param = process_options_structure( struct('sampling_method', 'empirical' ...
                                      ), param);

sampling_method = param.sampling_method;

model_mean = proteinModel.current_mean;

model_std = proteinModel.current_standard_deviation;

model_data = proteinModel.current_data;

% currently we assume there is no dependency between each voxels. 
switch sampling_method
    case 'gaussian'
        % sampling method Gaussian distribution.
        sampling_data = randn(size(model_mean)) .* model_std + model_mean;

    case 'empirical'
        % sampling method Empirical mean
        [f_ecdf, x_ecdf] = arrayfun(@(x) ecdf(model_data(:, x)), 1 : size(model_data, 2), 'UniformOutput', false);
        sampling_data = arrayfun(@(i) x_ecdf{i}(find(f_ecdf{i} < rand(), 1, 'last')), 1 : size(model_data, 2));        
end        

% normalize the image. 
sampling_data = sampling_data - min(sampling_data);

sampling_data = sampling_data ./ max(sampling_data(:));

template_image = proteinModel.t_cell_info.template_info.template_image;

I_raw = double(template_image); 

I_seg = template_image;

I_raw(I_raw > 0) = sampling_data;

end