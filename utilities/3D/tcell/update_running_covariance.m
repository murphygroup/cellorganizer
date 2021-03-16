function [result_sum_squares] = update_running_standard_deviation(given_preceding_mean, given_current_mean, given_preceding_sum_squares, given_parameters, given_index)
  % Update given_preceding_mean computed from (given_index - 1) data with given_parameters. Initial call should be with given_preceding_mean = zeros, given_index = 1, given_preceding_sum_squares = zeros.
  %
  % References:
  % Knuth, "The Art of Computer Programming, Volume 2: Seminumerical Algorithms", section 4.2.2
  % Welford 1962, "Note on a Method for Calculating Corrected Sums of Squares and Products"
  %
  % Tests:
  % d = 5;
  % a = randn(50, d) * randn(d, d);
  % am = mean(a, 1)
  % as = std(a, [], 1)
  % ac = cov(a)
  % arm = zeros(1, d);
  % ars = 0;
  % arc = zeros(d, d);
  % for index = 1:size(a, 1)
  %   pm = arm;
  %   cm = update_running_mean(pm, a(index, :), index);
  %   ps = ars;
  %   cs = update_running_standard_deviation(pm, cm, ps, a(index, :), index);
  %   pc = arc;
  %   cc = update_running_covariance(pm, cm, pc, a(index, :), index);
  %   arm = cm;
  %   ars = cs;
  %   arc = cc;
  % end
  % ars = finalize_running_standard_deviation(ars, size(a, 1));
  % arc = finalize_running_covariance(arc, size(a, 1));
  % arm
  % ars
  % arc
  % mean_error = norm(am - arm)
  % mean_standard_deviation_error = norm(as - ars)
  % mean_covariance_error = norm(ac - arc)
  %
  % Dependencies:
  % File Exchange:
  % Mine:
  % 
  % 2013-08-16 tebuck: Created.
  
  % Based on Welford 1962, "Note on a Method for Calculating Corrected Sums of Squares and Products" but in the form of Knuth's 4.2.2 (15, 16) used in update_running_standard_deviation:
  % result_sum_squares = given_preceding_sum_squares;
  % for row = 1:length(given_preceding_mean)
    % for column = row:length(given_preceding_mean)
      % % Direct Welford formula:
      % result_sum_squares(row, column) = given_preceding_sum_squares(row, column) + (given_parameters(row) - given_preceding_mean(row)) .* (given_parameters(column) - given_preceding_mean(column)) .* ((given_index - 1) / given_index);
      % % % Knuth-like:
      % % result_sum_squares(row, column) = given_preceding_sum_squares(row, column) + (given_parameters(row) - given_preceding_mean(row)) .* (given_parameters(column) - given_current_mean(column));
      % result_sum_squares(column, row) = result_sum_squares(row, column);
    % end
  % end
  % Direct Welford formula with vectorization:
  % result_sum_squares = given_preceding_sum_squares + (given_parameters - given_preceding_mean)' * (given_parameters - given_preceding_mean) .* ((given_index - 1) / given_index);
  % centralized_given_parameters = given_parameters - given_preceding_mean;
  % result_sum_squares = given_preceding_sum_squares + (centralized_given_parameters' * centralized_given_parameters) .* ((given_index - 1) / given_index);
  centralized_given_parameters_scaled = (given_parameters - given_preceding_mean) * sqrt((given_index - 1) / given_index);
  result_sum_squares = given_preceding_sum_squares + centralized_given_parameters_scaled' * centralized_given_parameters_scaled;
  % % Direct Welford formula with vectorization for multiple parameters:
  % centralized_given_parameters = given_parameters - repmat(given_preceding_mean, size(given_parameters, 1), 1)
  % result_sum_squares = given_preceding_sum_squares + ()' * spdiags(((given_index - 1:) / given_index), 0, length(current_standard_deviations), length(current_standard_deviations)) * (given_parameters - given_preceding_mean);
  
end


