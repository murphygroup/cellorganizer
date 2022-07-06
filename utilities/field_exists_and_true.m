function answer = field_exists_and_true( options, field, alltrue )
% FIELD_EXISTS_AND_TRUE Returns true if field exists and evaluates to true.
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
if alltrue
    answer = answer && options.(field);
else
    answer = answer && any(options.(field));
end

end
