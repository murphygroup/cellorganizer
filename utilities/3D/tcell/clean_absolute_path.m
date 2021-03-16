  function [result_path] = clean_absolute_path(given_path)
    % Remove any initial slash and reduce repeated slashes:
    result_path = given_path;
    % result_path = regexprep([result_path, '/'], '/+', '/');
    result_path = regexprep(strcat(result_path, '/'), '/+', '/');
  end
