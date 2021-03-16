function [protlabels,predlabels2,u,newallweights] = ...
    modeAcrossProtein( allweights, classlabels, antibodyids)

[vals, predlabels] = max(allweights,[],2);
u = unique(antibodyids);
protlabels = zeros(size(u));
predlabels2 = zeros(size(u));
tmp_idx = zeros(size(u));
newallweights = zeros(length(u),size(allweights,2));
for i=1:length(u)
    idx = find(antibodyids==u(i));
    newallweights(i,:) = sum(allweights(idx,:),1);
    protlabels(i) = classlabels(idx(1));
    ind = predlabels(idx);
    [val counts] = mode(ind);
    predlabels2(i) = val;
    if counts==length(ind)
        tmp_idx(i) = 1;
    end
end
idx = find(tmp_idx==0);
predlabels2(idx) = [];
protlabels(idx) = [];
u(idx) = [];
newallweights(idx,:) = [];
