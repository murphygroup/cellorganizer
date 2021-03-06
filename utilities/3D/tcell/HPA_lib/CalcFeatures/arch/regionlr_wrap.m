clear

addpath lib/classification
addpath lib/classification/logisticRegression
addpath lib/results
addpath lib/wrappers

load features/regionfeatures_all.mat

s = size(classlabels);
if s(1)==1
    classlabels = classlabels';
end

s = size(imagelist);
if s(1)==1
    imagelist = imagelist';
end

filterfeats

[r c] = find(isnan(allfeatures));
allfeatures(r,:) = [];
classlabels(r) = [];
antibodyids(r) = [];
imagelist(r) = [];

allfeatures = allfeatures(:,end-119:end);

[uu ii imagelabels] = unique(imagelist);
save( 'regionfeatures119.mat','allfeatures','classlabels','classes','antibodyids','imagelabels');

break
N = 5;
% determine CV splits
[proteins Ipr protlabels] = unique(antibodyids);

[proteins Ipr protlabels] = unique(antibodyids);
splits = partition( classlabels(Ipr), N, 13);
[trainidx, testidx] = partedsets( splits);

U = unique(classlabels);
s = zeros(size(U));

lambda = .3;
eta = 1e-3;
maxiter = 1000;
predlabels = zeros(length(classlabels),1);
allweights = zeros(length(classlabels),length(U));
for i=1:N
  disp(['evaluating fold: ' num2str(i)]);
    trainind = find(ismember(protlabels,trainidx{i}));
    testind = find(ismember(protlabels,testidx{i}));

    traindata = allfeatures( trainind,:);
    testdata = allfeatures( testind,:);
    
    trainlabels = classlabels(trainind);
    testlabels = classlabels(testind);
    
    % zscore standardization
    mu = mean(traindata,1);  st = std(traindata,[],1);
    traindata = (traindata - repmat( mu, [size(traindata,1) 1]))./ repmat( st, [size(traindata,1) 1]);
    testdata = (testdata - repmat( mu, [size(testdata,1) 1]))./ repmat( st, [size(testdata,1) 1]);

    % normalizing
    traindata = traindata./repmat( sqrt(sum(traindata(:,:).^2,2)), [1 size(traindata,2)]);
    testdata = testdata./repmat( sqrt(sum(testdata(:,:).^2,2)), [1 size(testdata,2)]);

    feat = [];
    for j=1:length(U)
        idx = find(trainlabels==U(j));
        feat{j} = traindata(idx,:);
        s(j) = length(idx);
    end
    idx_sda = ml_stepdisc( feat, ['sdaregionlr_' num2str(i) '.log']);
    idx = idx_sda(1:floor(length(idx_sda)/2));
    % this is overwritten below

if 1==2
    ma = 0;
    for j=1:1:length(idx_sda)
        idx_j = idx_sda(1:j);
        [lambda_best, eta_best, maxiter_best, ma_j] = parameterEstimationLR( trainlabels, traindata(:,idx_j), lambda, 1e-2, 100);

        if ma_j>ma
            ma = ma_j;
            idx = idx_j;
        end
    end
end

    w = lrtrain( traindata(:,idx), trainlabels, eta, lambda, maxiter);
    [predict_labels, probs] = lrtest( testdata(:,idx), w);

    [a ind] = max(probs,[],2);
    predlabels(testind) = predict_labels;
    allweights(testind,:) = probs;
end

[cc uu dd ww] = conmatrix( classlabels, predlabels);

disp(num2str(round(1000*uu)/10));
disp(' ');
disp(num2str(round(1000*dd)/10));

save('regionlrclass.mat', 'classlabels','allweights', 'classes', 'imagelist', 'antibodyids');
