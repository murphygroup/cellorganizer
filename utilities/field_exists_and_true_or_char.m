function answer = field_exists_and_true_or_char( options, field, alltrue )
% FIELD_EXISTS_AND_TRUE_OR_CHAR Returns true if field exists and either evaluates to true or is a char array.
%
% List Of Input Arguments  Descriptions
% -----------------------  ------------
% options                  A structure
% field                    A string specifying the field to check
% alltrue                  (optional) Boolean flag that indicates whether a true result requires the field to be all true (alltrue == true, Matlab's default behavior) or any true (alltrue == false). Default is true.

if nargin < 3
    alltrue = true;
end

answer = isfield(options,field);
answer = answer && true_or_char(options.(field), alltrue);
end
