function [t_cell_info, options] = tcell_adapt_models(t_cell_info, options)
% The function is called in img2model.m after building the model. It uses 
% t_cell_info and options as input and output. 
% 
% The purpose of the function is to adapt the model from the t cell pipeline 
% to the valid protein model in CellOrganizer. It setup the necessary options 
% for the model and put the learned model in the model structure. 
% 
% After this step, the model is a valid protein model, and will be further 
% process by functions the same as other protein models in CellOrganizer, 
% i.e. img2model.m, img2slml.m   
% 
% 2016-02-29: xruan: Copied from master_script_build_models.m.
% 
% 
% 
% 
% 
% % Dependencies:
% % From the File Exchange: export_fig, pmkmp, toolbox_graph, toolbox_fast_marching, affine, Snake3D, Mesh_voxelisation, dirr, inhull
% % From tebuck: modified version of Jieyue Li's HPA_lib, rasterize_mesh
%
% Author: Xiongtao Ruan
  
  
% Variables referenced more than once can be copied to local variables:

% get model type representation. 

t_cell_info = options.t_cell_info; 
[t_cell_info] = tcell_get_model_type_info(t_cell_info);

master_script_options = t_cell_info.options;

timepoints_to_include = options.timepoints_to_include;

model_types = t_cell_info.model_type_info.model_types;
number_model_types = t_cell_info.model_type_info.number_model_types;

% condition_sensor_combinations, number_condition_sensors, conditions, sensors, number_conditions, number_sensors
% model_types, number_model_types
% number_all_relative_times

if master_script_options.skip_model_building

warning('>>>> HACK, options.skip_model_building true!')

else

  % relevant_runs, beep, keyboard
  timepoints_to_include = t_cell_info.synapse_info.included_timepoints;
  for all_relative_time_index = 1 : numel(timepoints_to_include)
    relative_time = timepoints_to_include(all_relative_time_index);

    fprintf(' all_relative_time_index %d, relative_time % d\n', all_relative_time_index, relative_time)

    for model_type_index = 1:number_model_types
      % This is the .mat file, not the figure:
      current_model_filename = sprintf('%sreltime_%d.mat', master_script_options.model_prefix, relative_time);
        
      a = load([options.results_location, '/', current_model_filename]);
      % fixed bug in the case of process model multiple times
      if any(strcmp(fieldnames(a), 'model')) && numel(fieldnames(a)) <= 1
          disp('The model has been adpated, skip it!')
          continue;
      end
      model = struct();
      proteinModel = a;
      for fn = fieldnames(options.model)'
        proteinModel.(fn{1}) = options.model.(fn{1});
        model.(fn{1}) = options.model.(fn{1});
      end
      proteinModel.class = options.protein.class;
      proteinModel.type = options.protein.type;
      proteinModel.time_point = relative_time;
      proteinModel.name = master_script_options.model_prefix;
      
      model.proteinModel = proteinModel;

      save([options.results_location, '/', current_model_filename], 'model');             

    end 
  end
end

options.t_cell_info = t_cell_info;
  
end  

