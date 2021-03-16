function CompareTcellModels(model1, model2)
% compare the intensities of each voxel for two different models. 
% model1, and model2 are paths for the two models. 
% The figure will show the 3D slice of model1, model2, and the differences.
%
% Author: Xiongtao Ruan
% Date: Sep. 2016


debug = ~true;

if nargin < 2
    if debug
        model1 = 'LAT_reltime_0.mat';
        model2 = 'LAT_reltime_1.mat';        
    else
        error('Please provide two models for comparison!');
    end
end
try 
    a = load(model1);
    b = load(model2);
catch
    error('Unable to load model file!');
end

t_cell_info = a.model.proteinModel.t_cell_info;

template_crop_function = t_cell_info.template_info.template_crop_function;
model_reconstruction_functions = t_cell_info.model_type_info.model_reconstruction_functions; 
model_reconstruction_function = model_reconstruction_functions{1};

a_model_map = template_crop_function(model_reconstruction_function(a.model.proteinModel.current_mean));
b_model_map = template_crop_function(model_reconstruction_function(b.model.proteinModel.current_mean));

ab_model_diff_map = a_model_map - b_model_map;

close all;
figure,
set(gcf, 'Visible', 'on'), clf, 
figure_height = 250;
figure_width = 1500;
set(gcf, 'Position', [1, 1, figure_width, figure_height]);
set(gcf, 'Color', 'w')
image_to_show = [];
image_to_show = [image_to_show; a_model_map];
image_to_show = [image_to_show; b_model_map];
image_to_show = [image_to_show; ab_model_diff_map];

imagesc(reshape_2d(image_to_show));
axis off;
colorbar;

end