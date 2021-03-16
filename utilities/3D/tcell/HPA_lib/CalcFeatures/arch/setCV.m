function setCV( dtype, ftype, Nfolds)

if ~exist('dtype','var')
    dtype = 'field';
end
if ~exist('ftype','var')
    ftype = '9class';
    % ftype = '10class';
end

addpath lib/features
addpath lib/classification
addpath lib/classification/libsvm-mat-2.84-1
addpath lib/results
addpath lib/wrappers

load(['data/tmp/' dtype 'features_' ftype '.mat']);

if 1==2
idx = find(subfilters.cellabels==1);
subfeatures = subfeatures(idx,:);
sublabels.classlabels = sublabels.classlabels(idx);
sublabels.imagelist = sublabels.imagelist(idx);
sublabels.antibodyids = sublabels.antibodyids(idx);
end

if 1==1
    allfeatures = subfeatures;
    classlabels = sublabels.classlabels;
    classes = sublabels.classes;
    antibodyids = sublabels.antibodyids;
    imagelist = sublabels.imagelist;
else
    classes = sublabels.classes;
    [allfeatures,classlabels,antibodyids] = meddata( subfeatures, sublabels);
end

% allfeatures = allfeatures(:,[1:end-27 end-9]);
% allfeatures(:,[1:end-27 end-9]) = [];

if size(classlabels,1)==1
    classlabels = classlabels';
end

% determine CV splits
[proteins Ipr protlabels] = unique(antibodyids);
classlabels2 = zeros(size(classlabels));
for i=1:length(proteins)
    idx = find(proteins(i)==antibodyids);
    m = mode(classlabels(idx));
    classlabels2(idx) = m;
end
[c b] = hist(classlabels2(Ipr),1:1:max(classlabels2))
N = min(min(c),10);
if ~exist('Nfolds','var')
    Nfolds = 1:N;
end


% splits = partition( classlabels2(Ipr), N, 13);
splits = partition( classlabels, N, 13);
[trainidx, testidx] = partedsets( splits);

U = unique(classlabels);

c = 2.^(0:1:14);
g = [];
% g = 2.^(0:1:8);

predlabels = zeros(length(classlabels),1);
allweights = zeros(length(classlabels),length(U));
committee = [];
for i=Nfolds
  disp(['evaluating fold: ' num2str(i)]);
%     trainind = find(ismember(protlabels,trainidx{i}));
%     testind = find(ismember(protlabels,testidx{i}));
    trainind = trainidx{i};
    testind = testidx{i};

    traindata = allfeatures( trainind,:);
    testdata = allfeatures( testind,:);
    
    trainlabels = classlabels(trainind);
    testlabels = classlabels(testind);
    
    % zscore standardization
    [traindata, testdata, mu, st] = standardizeFeatures( traindata, testdata);

    % normalizing
    traindata = normalizeFeatures( traindata);
    testdata = normalizeFeatures( testdata);
    
    idx_feat = featureSelection( traindata, trainlabels, 'sda', num2str(i));
    idx = unique(idx_feat);
    wopt = classweights( trainlabels);

    if 1==2
        ma = 0;
        for j=21:10:length(idx_feat)
        disp(j);
            idx_j = idx_feat(1:j);
            [options_j,ma_j] = parameterEstimationSVM( trainlabels, traindata(:,idx_j), c, g, wopt);
        
            if ma_j>ma
                ma = ma_j;
                options = options_j;
                idx = idx_j;
            end
        end
    end
    
    % options = ['-s 0 -m 1000 -t 2 -c 64 -g 0.25 -b 1' wopt];
    options = ['-s 0 -m 1000 -t 0 -c 256 -b 1' wopt];
    % options = parameterEstimationSVM( trainlabels, traindata(:,idx), c, g, wopt);

    model = svmtrain( trainlabels, traindata(:,idx), options);
    [predict_label accuracy weights] = svmpredict( testlabels, testdata(:,idx), model, '-b 1');

    committee(i).model = model;
    committee(i).idx_feat = idx;
    committee(i).options = options;
    
    [a ind] = max(weights,[],2);
    predlabels(testind) = model.Label(ind);
    allweights(testind,model.Label) = weights;
end
[cc uu dd ww] = conmatrix( classlabels, predlabels)

Nstr = [];
for i=1:length(Nfolds)
    Nstr = [Nstr num2str(Nfolds(i)) '_'];
end
Nstr(end) = [];

save(['data/tmp/svm' dtype '_' ftype '_' Nstr '.mat'], 'allweights','committee');


proteins = unique(sublabels.antibodyids);
predlabels = zeros(size(proteins));
protlabels = zeros(size(proteins));
for i=1:length(proteins)
    idx = find(proteins(i)==sublabels.antibodyids);
    weights = sum(allweights(idx,:),1);
    [a label] = max(weights,[],2);
    
    predlabels(i) = label;
    
    protlabels(i) = mode(classlabels(idx));
end

[cc uu dd ww] = conmatrix( protlabels,predlabels)

