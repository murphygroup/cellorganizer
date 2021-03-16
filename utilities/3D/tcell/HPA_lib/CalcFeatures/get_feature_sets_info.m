function [all_feature_sets, feature_length_dictionary, required_channels_dictionary, original_feature_set_indices] = get_feature_sets_info()
  % Returns all feature set names, their feature counts (in a
  % dictionary formatted as a structure), the default sets as a
  % list of indices, and the channels needed to compute them.
  % 
  % Created 2011-12-11 tebuck.

original_feature_sets = {...
    'mutualInfo', ...
    'mutualInfox2', ...
    'nonObjFluor', ...
    'nonObjFluorx2', ...
    'obj', ...
    'objx2', ...
    'overlap', ...
    'overlapx2', ...
    'texture', ...
    'texturex2', ...
    'texturex4', ...
    'tas', ...
    'tasx2', ...
    'objRegion',...
    'objRegionx2',...
};
original_number_features = [30, 30, 4, 4, 18, 18, 39, 39, 156, 156, 156, 18, 18, 14, 14];
% tebuck: custom features:
all_feature_sets = original_feature_sets; 
all_number_features = original_number_features; 
all_feature_sets = [all_feature_sets, {'nucStats'}]; 
all_number_features = [all_number_features, 20]; 
all_feature_sets = [all_feature_sets, {'protTotalIntensity'}]; 
all_number_features = [all_number_features, 1]; 
original_feature_set_indices = 1:length(original_feature_sets); 
%required_channels = cell(length(original_feature_sets), 1);
required_channels = repmat({{}}, length(all_feature_sets), 1);
%required_channels([1, 2, 7, 8]) = repmat({{'protein', 'nuclear', 'tubulin', 'er'}}, 4, 1); 
required_channels([1, 2, 3, 4, 7, 8]) = {{'protein', 'nuclear', 'tubulin', 'er'}}; 
required_channels([5, 6, 14, 15]) = {{'protein', 'nuclear'}}; 
required_channels([9, 10, 11, 12, 13, 17]) = {{'protein'}}; 
required_channels([16]) = {{'nuclear'}}; 


% Build a dictionary of number of features per feature set:
all_number_features_cell = num2cell(all_number_features);
feature_length_dictionary = struct();
for field_index = 1:length(all_feature_sets)
  feature_length_dictionary.(all_feature_sets{field_index}) = all_number_features(field_index);
end
%feature_length_dictionary

required_channels_dictionary = struct();
for field_index = 1:length(all_feature_sets)
  required_channels_dictionary.(all_feature_sets{field_index}) = required_channels{field_index};
end
%required_channels_dictionary