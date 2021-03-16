function [options] = tcell_img2param(options)
% The script is called in tcell_imgs2param.m. The script calls tcell_to_parameter.m to 
% extract parameters for each cell.
% 
% Author: Xiongtao Ruan


cellfit.type = 'morphing';
cell_option = struct();
cell_option.cell_index = options.cell_index;
cell_param = tcell_to_parameter(cell_option, options);   
% because the need to save options, there is no need to
% save t_cell_info for every cell. 
rmfield(options, 't_cell_info');

end