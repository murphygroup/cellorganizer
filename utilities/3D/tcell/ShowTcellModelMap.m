function ShowTcellModelMap(model_filename, mode)
% The script shows slices through the 3D map for a specific time point. 
%
% Author: Xiongtao Ruan
% Date: Sep. 2016
% 10/05/2018 add option to chose whether show model mean or standard deviation. 


debug = ~true;
if nargin < 2
	mode = 'mean';   % 'mean' or 'std'
end

if nargin < 1
    if debug
        model_filename = 'LAT_reltime_1.mat';
    else
        error('Please provide model filename!');
    end
end
try 
    load(model_filename);
catch
    error('Unable to load model file!');
end

t_cell_info = model.proteinModel.t_cell_info;
template_crop_function = t_cell_info.template_info.template_crop_function;
model_reconstruction_functions = t_cell_info.model_type_info.model_reconstruction_functions; 
 
model_reconstruction_function = model_reconstruction_functions{1};

switch mode
case 'mean'
	current_model_values = model.proteinModel.current_mean;
case 'std'
	current_model_values = model.proteinModel.current_standard_deviation;
end

model_map = model_reconstruction_function(current_model_values);

figure,
set(gcf, 'Visible', 'on'), clf, 
figure_height = 100;
figure_width = 1500;
set(gcf, 'Position', [1, 1, figure_width, figure_height]);
set(gcf, 'Color', 'w')
imagesc(reshape_2d(template_crop_function(model_map)));
h = colorbar('FontSize',6);
h.Position = [0.92 0.2 0.02 0.65];
axis off;
% colorbar('FontSize',5, 'Ticks', [0,  max(model_map(:)) * 0.6]);

end
