clear

addpath lib/classification
addpath lib/classification/libsvm-mat-2.84-1
addpath lib/results

load features/regionfeatures_all.mat

filtertrainfeats_script


[r c] = find(isnan(allfeatures));
allfeatures(r,:) = [];
classlabels(r) = [];
antibodyids(r) = [];
imagelist(r) = [];



N = 5;
[proteins Ipr protlabels] = unique(antibodyids);
splits = partition( classlabels(Ipr)', N, 13);
[trainidx, testidx] = partedsets( splits);


U = unique(classlabels);
s = zeros(size(U));

c = 2.^[0:2:8];
g = 2.^[-3:1:2];
c = c(4);
g = g(2);

predlabels = zeros(length(classlabels),1);
allweights = zeros(length(classlabels),max(classlabels));
for i=1:N
tic;
    trainind = find(ismember(protlabels,trainidx{i}));
    testind = find(ismember(protlabels,testidx{i}));
    
    traindata = allfeatures(trainind,:);
    testdata = allfeatures(testind,:);
    trainlabels = classlabels(trainind);
    testlabels = classlabels(testind);
    
    % zscore standardization
    mu = mean(traindata,1);  st = std(traindata,[],1);
    traindata = (traindata - repmat( mu, [size(traindata,1) 1]))./ repmat( st, [size(traindata,1) 1]);
    testdata = (testdata - repmat( mu, [size(testdata,1) 1]))./ repmat( st, [size(testdata,1) 1]);

    % normalizing
    traindata = traindata./repmat( sqrt(sum(traindata(:,:).^2,2)), [1 size(traindata,2)]);
    testdata = testdata./repmat( sqrt(sum(testdata(:,:).^2,2)), [1 size(testdata,2)]);

    % [traindata,testdata] = featnorm( traindata, testdata);
    % traindata = 2*traindata-1;
    % testdata = 2*testdata-1;
    
    feat = [];
    for j=1:length(U)
         feat{j} = traindata(find(trainlabels==U(j)),:);
    end
    idx_sda = ml_stepdisc( feat,['sdaregionsvm_' num2str(i) '.log']);
    idx_sda = idx_sda(1:60); % idx_sda(1:floor(length(idx_sda)/2));

    wopt = classweights( trainlabels);
    options = ['-s 0 -m 1000 -t 2 -c 64 -g 0.25 -b 1' wopt];

    % [model,options] = trainClassifier( trainlabels, traindata(:,idx_sda), c, g);
    model = svmtrain( trainlabels', traindata(:,idx_sda), options);
    [predict_label accuracy weights] = svmpredict( testlabels', testdata(:,idx_sda), model, '-b 1');
    
    [a ind] = max(weights,[],2);
    predlabels(testind) = model.Label(ind);
    allweights(testind,model.Label) = weights;
toc
end

save('regionclass.mat','classlabels','allweights','predlabels');
