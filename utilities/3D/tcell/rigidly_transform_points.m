function [transformed_points, transform, transform_affine] = rigidly_transform_points(points, given_volume_size, given_rotation, given_translation, given_rotation_center, use_inverse_rotation)
  % Rigidly rotate and then translate a set of 3D points. given_rotation is [azimuth, elevation] where azimuth is in the XY plane. given_translation is [x, y, z] where x is horizontal/column. Rotation is currently about given_volume_size / 2 to be consistent with rigidly_transform_image by default, but this point can be modified with argument given_rotation_center. For convenience, given_volume_size is YXZ like returned by size(volume_image). If use_inverse_rotation is true, the transformation implied by given_rotation is inverted when applied to undo a spherical coordinates rotation.
  % 
  % 2013-04-07 tebuck: copied from align_time_series.m.
  % 2013-04-28 tebuck: modified to have optional parameter given_rotation_center to use instead of the center of the volume.
  % 
  % Dependencies:
  % From the File Exchange: affine
  % From tebuck: rotation_matrix3h, translation_matrix3h
  
  if ~exist('given_rotation_center', 'var') || isempty(given_rotation_center)
    % Image size is used for centering by default:
    given_rotation_center = given_volume_size / 2 + .5;
    given_rotation_center = given_rotation_center([2, 1, 3]);
  end
    
  if ~exist('use_inverse_rotation', 'var')
    use_inverse_rotation = false;
  end
    
  translation_component = given_translation;
  rotation_component = rotation_matrix3h(3, given_rotation(1)) * rotation_matrix3h(2, -given_rotation(2));
  if use_inverse_rotation
    rotation_component = pinv(rotation_component);
  end
  
  % Set up total transform:
  transform = translation_matrix3h(translation_component) * translation_matrix3h(given_rotation_center) * rotation_component * translation_matrix3h(-(given_rotation_center));
  transform_affine = translation_matrix3h((translation_component - 1)) * translation_matrix3h((given_rotation_center - 1)) * rotation_component * translation_matrix3h(-(given_rotation_center - 1));
  
  % Transform current_landmarks in homogeneous coordinates:
  transformed_points = [points, ones(size(points, 1), 1)];
  % whos transform transformed_points
  transformed_points = transform * transformed_points';
  transformed_points = transformed_points(1:3, :)';
end
