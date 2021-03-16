function ts = tz_expandcombabbrv(funname)
%TZ_EXPANDCOMBABBRV Expand combination of abbreviations
%   TS = TZ_EXPANDCOMBABBRV(S) returns full words of S, which is the
%   combination of a set abbreviations
%   It is useful for taking a quick look at how the function or a variable
%   is named :).
%   TS is cell array.

%   12-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argumentis required')
end

nameLength = length(funname);
loopTime = pow2(nameLength-1);
penalty = zeros(1,loopTime);

for i = 1:loopTime
    wordPattern = dec2bin(loopTime+i-1);
    abbreviationArray = tz_strtok_pos(funname,wordPattern);
    for j = 1:length(abbreviationArray)
        if strcmp(tz_expandabbrv( abbreviationArray{j} ), 'unknown word')
            penalty(i) = penalty(i) + length(abbreviationArray{j});
        end
    end
    penalty(i) = penalty(i) + log(length(abbreviationArray));
    if penalty(i) < 1
        penalty(end+1) = 0;
        penalty(i+1:end) = [];
        break;
    end
end

[leastPenaly,bestPattern] = min(penalty);
bestPattern = dec2bin(loopTime+bestPattern-1);

ts = tz_strtok_pos(funname,bestPattern);

for i = 1:length(ts)
    fullWord = tz_expandabbrv(ts{i});
    if ~strcmp(fullWord, 'unknown word')
        ts{i} = fullWord;
    end
end



