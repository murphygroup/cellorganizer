function [trainidx,testidx, protlabels,N] = splitbyproteins( classlabels, abids2, Ninit)

if ~exist('Ninit','var')
    Ninit = 10;
end


% determine CV splits
[proteins Ipr protlabels] = unique(abids2);
classlabels2 = zeros(size(classlabels));
for i=1:length(proteins)
    idx = find(proteins(i)==abids2);
    m = mode(classlabels(idx));
    classlabels2(idx) = m;
end
%[c b] = hist(classlabels2(Ipr),1:1:max(classlabels2))
[c b] = hist(classlabels2(Ipr),unique(classlabels2))
% [c b] = hist(classlabels,1:1:max(classlabels))



N = min(min(c),Ninit);

%keyboard

splits = partition( classlabels2(Ipr), N, 13);
% splits = partition( classlabels, N, 13);
[trainidx, testidx] = partedsets( splits);
