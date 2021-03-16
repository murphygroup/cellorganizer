function [trainidx,testidx] = partedsets( splits)
% [TRAINIDX, TESTIDX] = PARTEDSETS( SPLITS)
%   index data splits into N training and testing sets

trainidx = [];
testidx = [];
trainidx{length(splits)} = [];
testidx{length(splits)} = [];
for i=1:length(splits)
    idx = 1:length(splits);
    idx(i) = [];
    testidx{i} = splits{i};
    for j=1:length(idx)
        trainidx{i} = [trainidx{i}; splits{idx(j)}];
    end
end

