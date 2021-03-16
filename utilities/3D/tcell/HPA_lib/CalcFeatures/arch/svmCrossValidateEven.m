function committee = svmCrossValidateEven( allfeatures, classlabels, antibodyids, randS, dtype, rratio) 
% function [protlabels,predlabels,proteins,newallweights] = ...
%     svmCrossValidate( allfeatures, classlabels, antibodyids, randS)

addpath lib/classification
addpath lib/classification/libsvm-mat-2.84-1
addpath lib/results





% determine CV splits
[proteins Ipr protlabels] = unique(antibodyids);
classlabels2 = zeros(size(classlabels));
for i=1:length(proteins)
    idx = find(proteins(i)==antibodyids);
    m = mode(classlabels(idx));
    classlabels2(idx) = m;
end
[c b] = hist(classlabels2(Ipr),1:1:max(classlabels2));
N = min(min(c),10);
if ~exist('Nfolds','var')
    Nfolds = 1:N;
end

rand('seed',randS);

% mi = min(c);
mi = floor(min(c)*rratio);
mii = [];
for i=1:length(c)
    rp = randperm(c(i));
    idx = find(classlabels2(Ipr)==i);
    mii = [mii; idx(rp(1:mi))];
end

protlist = proteins(mii);

idx = [];
for i=1:length(protlist)
    idx = [idx; find(protlist(i)==antibodyids)];
end

allfeatures = allfeatures(idx,:);
classlabels = classlabels(idx);
antibodyids = antibodyids(idx);







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

splits = partition( classlabels2(Ipr), N, 13);

U = unique(classlabels);



% U = unique(classlabels);
% [c b] = hist(classlabels,1:1:max(U));
% N = min(10,min(c));
% splits = partition( classlabels, N, 13);


[trainidx, testidx] = partedsets( splits);
predlabels = zeros(length(classlabels),1);
allweights = zeros(length(classlabels),length(U));


c = 2.^(0:1:12);
g = 2.^(0:1:8);
committee = [];
for i=1:N
    trainind = find(ismember(protlabels,trainidx{i}));
    testind = find(ismember(protlabels,testidx{i}));
%     trainind = trainidx{i};
%     testind = testidx{i};

    traindata = allfeatures( trainind,:);
    testdata = allfeatures( testind,:);
    
    trainlabels = classlabels(trainind);
    testlabels = classlabels(testind);
    
    % zscore standardization
    [traindata, testdata, mu, st] = standardizeFeatures( traindata, testdata);

    % normalizing
    traindata = normalizeFeatures( traindata);
    testdata = normalizeFeatures( testdata);
    
    ut = clock;
    ut = [num2str(ut(end-2)) num2str(ut(end-1)) num2str(ut(end))];
    ut(ut=='.') = '_';
    ut = [ut '_' num2str(randS)];

    idx_feat = featureSelection( traindata, trainlabels, 'sda', ut);

    wopt = classweights( trainlabels);
    options = ['-s 0 -m 1000 -t 2 -c 2048 -g 4 -b 1' wopt];
    options = parameterEstimationSVM( trainlabels, traindata(:,idx_feat), c, g, wopt);

    model = svmtrain( trainlabels, traindata(:,idx_feat), options);
%     [predict_label accuracy weights] = svmpredict( testlabels, testdata(:,idx_feat), model, '-b 1');

    committee(i).model = model;
    committee(i).idx_feat = idx_feat;
    committee(i).options = options;
    committee(i).mu = mu;
    committee(i).st = st;
    
%     [a ind] = max(weights,[],2);
%     predlabels(testind) = model.Label(ind);
%     allweights(testind,model.Label) = weights;
end

% [cc uu dd ww] = conmatrix( classlabels, predlabels)

% save(['data/tmp/svmEven_' dtype '_' num2str(randS) '.mat'],'committee','cc','uu','dd','ww');
