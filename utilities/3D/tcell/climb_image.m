function [positions] = climb_image(landscape, positions)
% 2013-01-21 tebuck: Extracted from
% segment_using_snakes.m.
% Dependencies:
% tebuck: smooth_gradient
  
  % Improve seed locations using nearby maxima:
  
  number_positions = size(positions, 1);

  [gx, gy, gz] = tcell_smooth_gradient(landscape, 'scharr5', false, true);
  
  % Simple local gradient descent:
  % step_size = 0.1;
  % step_size = 0.25;
  step_size = 0.5;
  number_previous_positions_to_use = 10;
  convergence_ratio = 0.1;
  % convergence_ratio = 0.3;
  old_displacements = zeros(number_positions, 3, 0);
  % This gets stuck, just terminate instead of RK integrator:
  % maximum_iterations = 100;
  maximum_iterations = 100;

  % old_positions = positions + 10.;
  old_positions = positions + convergence_ratio * step_size * 2;
  old_displacements = zeros(number_positions, 3, 0);
  
  % refinement_start_time = tic;
  running = true;
  while running
    old_positions(:, :, end + 1) = positions;
    x = positions(:, 1);
    y = positions(:, 2);
    z = positions(:, 3);
    % positions
    displacements = [interp3(gx, x, y, z), interp3(gy, x, y, z), interp3(gz, x, y, z)];
    displacement_sizes = repmat(sqrt(sum(displacements.^2, 2)), 1, 3);
    displacements = displacements ./ displacement_sizes * step_size;
    old_displacements(:, :, end + 1) = displacements;
    % if any(displacement_sizes > step_size)
      % displacements(displacement_sizes > step_size, :) = displacements(displacement_sizes > step_size, :) ./ displacement_sizes(displacement_sizes > step_size, :) * step_size;
    % end
    positions = positions + step_size * displacements;
    % Keep points inside the image:
    positions(positions < 1) = 1;
    positions(positions(:, 1) > size(gx, 2), 1) = size(gx, 2);
    positions(positions(:, 2) > size(gx, 1), 2) = size(gx, 1);
    positions(positions(:, 3) > size(gx, 3), 3) = size(gx, 3);
    % positions
    if size(old_positions, 3) > number_previous_positions_to_use
      % previous_positions = old_positions(:, :, end - number_previous_positions_to_use + 1:end);
      previous_displacements = old_displacements(:, :, end - number_previous_positions_to_use + 1:end);
      % running = mean(mean(mean(previous_displacements, 3) ./ std(previous_positions, 3))) > convergence_ratio;
      % running = mean(previous_displacements(:)) > convergence_ratio * step_size;
      running = max(previous_displacements(:)) > convergence_ratio * step_size;
      % % Use the variance of a uniform distribution as a convergence criterion for steps oscillating around a point:
      % % uniform_variance = (1/12 * (2 * step_size)^2)
      % uniform_variance = (1/12 * (step_size)^2);
      % mean_current_variance = var(previous_displacements, 0, 3);
      % mean_current_variance = mean(mean_current_variance(:));
      % % mean_current_variance = sqrt(sum(previous_displacements.^2, 2)), pause
      % % mean_current_variance = mean(mean_current_variance);
      % compared_values = [abs(mean_current_variance - uniform_variance), convergence_ratio * step_size]
      % % running = running && norm(cov(previous_displacements) - eye(3) * uniform_variance) > convergence_ratio * step_size;
      % running = running && compared_values(1) > compared_values(2);
    end
    running = running && (maximum_iterations > size(old_displacements, 3));
  end
  % iterations_run = size(old_displacements, 3)
  % refinement_time = toc(refinement_start_time)
  
  % imagesc(mean(landscape, 3)); colormap(gray); hold on; plot3(squeeze(old_positions(:, 1, :))', squeeze(old_positions(:, 2, :))', squeeze(old_positions(:, 3, :))', 'LineWidth', 2); hold off; pause


end