clear

setup
addpath lib/classification
addpath lib/classification/libsvm-mat-2.84-1


parse_wrap

% classlabels
% imagelist

ulist = listunixformat( imagelist);

featsets = {'overlapField','tasField','texture128Field'};
numfeats = [9, 18, 26];
encount = cumsum(numfeats);
stcount = [1 encount(1:end-1)+1];

for i=1:length(featsets)
    tmplist = findreplacestring( ulist, '/','_');
    tmplist = findreplacestring( tmplist, '_images_HPA_confocal_', 'data/features/');
    flists{i} = findreplacestring( tmplist, '.tif', ['_' featsets{i} '.mat']);
    featlists{i} = listmatlabformat( flists{i})';
end

allfeatures = zeros(length(imagelist),encount(end));
idx = zeros(length(imagelist),1);
for i=1:length(imagelist)
    breaker = 0;
    features = zeros(1,encount(end));
    for j=1:length(featsets)
        if ~exist(featlists{j}{i},'file')
            breaker = 1;
            continue;
        end
        load(featlists{j}{i});
        features(stcount(j):encount(j)) = feats;
    end
    
    if ~breaker
        idx(i) = 1;
        allfeatures(i,:) = features;
    end
end

ind = find(idx);
allfeatures = allfeatures(ind,:);
imagelist = imagelist(ind);
classlabels = classlabels(ind);
antibodyids = antibodyids(ind);

[classes Icl classlabels] = unique(classlabels);
[proteins Ipr protlabels] = unique(antibodyids);

[c b] = hist( classlabels(Ipr), [1:1:max(classlabels(:))]);

idx_rmclass = find(c<10);
idx = [];
for i=1:length(idx_rmclass)
    idx = [idx; find(idx_rmclass(i)==classlabels)];
end

classlabels(idx) = [];
imagelist(idx) = [];
allfeatures(idx,:) = [];
antibodyids(idx) = [];

[proteins Ipr protlabels] = unique(antibodyids);

