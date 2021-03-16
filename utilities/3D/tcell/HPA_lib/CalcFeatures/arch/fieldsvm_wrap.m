function fieldsvm_wrap(doN)
% clear

addpath lib/classification
addpath lib/classification/libsvm-mat-2.84-1
addpath lib/results
addpath lib/wrappers


savepathname = ['fieldclass_' num2str(doN) '.mat'];
if exist( savepathname, 'file')
    delete(savepathname);
end

load data/tmp/fieldfeatures_all.mat
filterfeats

N = 10;
% determine CV splits
[proteins Ipr protlabels] = unique(antibodyids);
% splits = partition( classlabels(Ipr), N, 13);
% [trainidx, testidx] = partedsets( splits);

splits = partition( classlabels, N, 13);
[trainidx, testidx] = partedsets( splits);


U = unique(classlabels);
s = zeros(size(U));

c = [2 4 8 16 32 64 128 1024];
g = [0.125 0.25 0.5 1];

c = 2.^(-2:1:12);
g = 2.^(-8:1:2);
% c = 64;
% g = 0.25;

predlabels = zeros(length(classlabels),1);
allweights = zeros(length(classlabels),length(U));
for i=doN %i=1:N
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
    idx_sda = ml_stepdisc( feat, ['sdafield_' num2str(i) '.log']);
    idx = idx_sda(1:floor(length(idx_sda)/2));
    % this is overwritten below

    wopt = classweights( trainlabels);
    
if 1==1
ma = 0;
    for j=1:length(idx_sda)
        idx_j = idx_sda(1:j);
        [options_j,ma_j] = parameterEstimation( trainlabels, traindata(:,idx_j), c, g, wopt);
        
        if ma_j>ma
            ma = ma_j;
            options = options_j;
            idx = idx_j;
        end
    end
end

    % options = ['-s 0 -m 1000 -t 2 -c 64 -g 0.25 -b 1' wopt];
    
    model = svmtrain( trainlabels, traindata(:,idx), options);
    [predict_label accuracy weights] = svmpredict( testlabels, testdata(:,idx), model, '-b 1');

    [a ind] = max(weights,[],2);
    predlabels(testind) = model.Label(ind);
    allweights(testind,model.Label) = weights;
end
[cc uu dd ww] = conmatrix( classlabels, predlabels)

% save('fieldclass.mat');
save(savepathname, 'classlabels','allweights','testind','model','weights');
