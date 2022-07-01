function answer = true_or_char( value, alltrue )
% TRUE_OR_CHAR Returns true if value either evaluates to true or is a char array.
%
% List Of Input Arguments  Descriptions
% -----------------------  ------------
% value                    A string specifying the field to check
% alltrue                  (optional) Boolean flag that indicates whether a true result requires the field to be all true (alltrue == true, Matlab's default behavior) or any true (alltrue == false). Default is true.

if nargin < 2
    alltrue = true;
end

if alltrue
    answer = ~isempty(value) && (ischar(value) || value);
else
    answer = ~isempty(value) && (ischar(value) || any(value));
end

end
