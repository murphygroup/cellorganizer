function str2 = tz_strconvert(str1,setname1,setname2)
%TZ_STRCONVERT Convert from one string to the other.
%   STR2 = TZ_STRCONVERT(STR1,SET1,SET2) returns a cell arry of strings,
%   which are strings in a set with the name SET2 corresponding to the
%   the strings in SET1. STR1 must be a string in SET1, otherwise an
%   empty cell array will be returned. Currectlly the available string
%   sets are:
%       location - 2dhela location names
%       pattern - 2dhela pattern names
%   
%   See also

%   25-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

locationSet = {'DNA','ER','Golgi','Golgi','lysosome','mitochondria', ...
        'nucleolus','transferrin receptor','actin','tubulin'};
patternSet = {'DAPI','ER','giant','gpp130','LAMP2','mito.','nucleolar', ...
        'TfR','phal','tubul'};

inputSetNames = {setname1,setname2};

for i=1:2
    switch inputSetNames{i}
    case 'location'
        inputSets{i} = locationSet;
    case 'pattern'
        inputSets{i} = patternSet;
    end
end

index1 = strmatch(str1,inputSets{1});
if ~isempty(index1)
    str2 = inputSets{2}(index1);
else
    str2 = [];
end