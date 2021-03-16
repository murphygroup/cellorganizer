function [protlabels,predlabels,proteins,newallweights] = ...
    voteAcrossProtein( allweights, classlabels, antibodyids)

proteins = unique(antibodyids);
predlabels = zeros(size(proteins));
protlabels = zeros(size(proteins));
newallwegihts = zeros(length(proteins),size(allweights,2));
for i=1:length(proteins)
    idx = find(proteins(i)==antibodyids);
    weights = sum(allweights(idx,:),1);
    [a label] = max(weights,[],2);
    
    predlabels(i) = label;
    
    protlabels(i) = mode(classlabels(idx));
    newallweights(i,:) = weights/length(idx);
end

[cc uu dd ww] = conmatrix( protlabels,predlabels);
