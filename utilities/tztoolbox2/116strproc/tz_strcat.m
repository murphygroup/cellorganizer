function ss = tz_strcat(s1,ss2)
%TZ_STRCAT Concatenate strings
%   SS = TZ_STRCAT(S1,SS2) horizontally concatenates the string S1 and the
%   [string array] SS2. The returned value is a [string array] with each
%   element to be the concatenated string of S1 and the correpsonding
%   element in SS2.
%   
%   See also

%   29-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

for i=1:length(ss2)
    ss{i} = strcat(s1,ss2{i});
end

if size(ss2,1)>1
    ss = ss';
end