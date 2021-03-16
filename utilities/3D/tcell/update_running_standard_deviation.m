function [result_standard_deviation] = update_running_standard_deviation(given_preceding_mean, given_current_mean, given_preceding_standard_deviation, given_parameters, given_index)
  % Update given_preceding_mean computed from (given_index - 1) data with given_parameters. Initial call should be with given_preceding_mean = zeros, given_index = 1, given_preceding_standard_deviation = zeros.
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
  
  result_standard_deviation = given_preceding_standard_deviation + (given_parameters - given_preceding_mean) .* (given_parameters - given_current_mean);
  
end


