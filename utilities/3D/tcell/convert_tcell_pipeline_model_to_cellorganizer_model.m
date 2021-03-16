function convert_tcell_pipeline_model_to_cellorganizer_model(model_name, output_dir, reference_model)
% convert models from the t cell pipeline to valid cellorganizer model
% 
% Author: Xiongtao Ruan
% Date: Aug. 31, 2017

if nargin < 3 
    reference_model = './models_1/LAT_reltime_-2.mat';
end

if nargin < 1
    model_name = 'Full Stimulus_LAT_reltime  0_standardized_voxels.mat';
end

if nargin < 2
    output_dir = pwd;
end

a = load(reference_model);

b = load(model_name);
pattern = 'reltime +(-?\d+)_';
reltime = regexp(model_name, pattern, 'tokens');

if ~isempty(reltime)
    reltime = str2double(reltime{1}{1});
end

[pathstr, model_file_name] = fileparts(model_name);

prefix_pattern = '^([^_]+_[^_]+_)reltime';
prefix_cell = regexp(model_file_name, prefix_pattern, 'tokens');

if ~isempty(prefix_cell)
    prefix = prefix_cell{1}{1};
end

model = a.model;

output_name = sprintf('%sreltime_%d.mat', prefix, reltime);

model.proteinModel.name = prefix;
model.proteinModel.time_point = reltime;
model.proteinModel.relative_time = reltime;

model.proteinModel.current_mean = b.current_mean;
model.proteinModel.current_standard_deviation = b.current_standard_deviation;
model.proteinModel.current_data = b.current_data;

save([output_dir, '/', output_name], 'model');


end