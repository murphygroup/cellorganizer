function [t_cell_info, options] = tcell_imgs2params(t_cell_info, options)
% The script is called in img2model.m.
% the purpose of the function is to s to first parameterize the images before
% building the model. This script is the main script for image parameterization.

synapse_info = t_cell_info.synapse_info;
all_run_data = synapse_info.all_run_data;
num_imgs = size(all_run_data, 1);

% looping for each image for parameterization
for i = 1:num_imgs
    paramfiles{i} = [options.paramdir filesep 'param' num2str(i) '.mat'];
    
    tmpfile = [paramfiles{i} '.tmp'];
    if exist(tmpfile, 'file')
        continue
    elseif exist(paramfiles{i}, 'file')
        isdone(i) = true;
        continue
    end
    
    system(['touch "' tmpfile '"']);
    
    
    savedir = [options.paramdir filesep 'param' num2str(i)];
    
    options.cell_index = i;
    
    % call the function to extract parameter.
    [cell_params] = tcell_img2param(options);
    save(paramfiles{i}, '-struct', 'cell_params')
    isdone(i) = true;
    
    delete(tmpfile);
    try
    catch the_error
        warning('Unable to extract parameters.')
        disp( 'Check the images exist or that you are using the correct options.' );
    end
end
end