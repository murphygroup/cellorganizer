function [transformed_image] = rigidly_transform_image(given_image, given_rotation, given_translation)
  % Rigidly rotate and then translate a 3D image. given_rotation is [phi, psi], given_translation is [x, y, z] where x is horizontal/column. Rotation is currently about size(given_image) / 2.
  %
  % 2013-04-07 tebuck: copied from align_time_series.m.
  % Dependencies:
  % From the File Exchange: affine
  % From tebuck: rotation_matrix3h, translation_matrix3h
  
  % pivot_component = size(given_image) / 2;
  pivot_component = size(given_image) / 2 + .5;
  pivot_component = pivot_component([2, 1, 3]);
  
  % Cumulative transformations:
  rotation_component = rotation_matrix3h(3, given_rotation(1)) * rotation_matrix3h(2, -given_rotation(2));

  % % Use rotation vectors with length indicating number of revolutions about the axis which is the direction:
  % rotation_vector = given_rotation;
  % rotation_revolutions = norm(rotation_vector);
  % if rotation_revolutions > 0
    % rotation_vector = rotation_vector ./ rotation_revolutions;
  % else
    % rotation_vector = [1, 0, 0];
  % end
  % rotation_component = rotation_matrix3h(rotation_vector, rotation_revolutions * 360);
  
  % Set up transform for current_image:
  transform_affine = translation_matrix3h(given_translation - 1) * translation_matrix3h(pivot_component - 1) * rotation_component * translation_matrix3h(-(pivot_component - 1));
  
  % Use YXZ order for affine:
  transform_affine = transform_affine([2, 1, 3, 4], [2, 1, 3, 4]);
  transformed_image = affine(given_image, transform_affine, [], false, [], [], 'crop');
  % transformed_image = permute([2, 1, 3]);
end

