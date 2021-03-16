clear

addpath lib/classification
addpath lib/classification/libsvm-mat-2.84-1
addpath lib/results
addpath lib/wrappers

load features/fieldfeatures_all.mat

filterfeats

N = 10;
% determine CV splits
[proteins Ipr protlabels] = unique(antibodyids);
splits = partition( classlabels(Ipr), N);
[trainidx, testidx] = partedsets( splits);

U = unique(classlabels);
s = zeros(size(U));
numk = round(median(s));
numk = 5;

predlabels = zeros(length(classlabels),1);
distances = zeros(length(classlabels),numk);
indices = zeros(length(classlabels),numk);
for i=1:N
    trainind = find(ismember(protlabels,trainidx{i}));
    testind = find(ismember(protlabels,testidx{i}));

    traindata = allfeatures( trainind,:);
    testdata = allfeatures( testind,:);
    
    trainlabels = classlabels(trainind);
    
    
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
    idx = ml_stepdisc( feat);

    subN = 3;
    subnumk = numk;
    subsplits = partition( trainlabels, subN);
    [subtrainidx, subtestidx] = partedsets( subsplits);
    
    val = zeros(length(idx),1);
    for k=1:1:length(idx)
    disp(k);
        idx_sda = idx(1:k);
        subpredlabels = zeros(length(trainlabels),subnumk);
        % subdistances = zeros(length(trainlabels),subnumk);
        % subindices = zeros(length(trainlabels),subnumk);
        for j=1:subN
            subtraindata = traindata( subtrainidx{j},:);
            subtestdata = traindata( subtestidx{j},:);
    
            subtrainlabels = trainlabels( subtrainidx{j});
        
            knn = pdist([subtestdata(:,idx_sda); subtraindata(:,idx_sda)], 'euclidean');
            knn = squareform(knn);
            knni = knn(1:size(subtestdata(:,idx_sda),1),size(subtestdata(:,idx_sda),1)+1:end);

            [y ind] = sort(knni,2);
            nlabels = subtrainlabels(ind(:,1:subnumk));
            % subindices(subtestidx{j},:) = nlabels;
            % subdistances(subtestidx{j},:) = y(:,1:subnumk);
            subpredlabels(subtestidx{j}) = mode(nlabels,2);
        end
        [cc uu dd ww] = conmatrix( trainlabels, subpredlabels);
        disp([ uu ww]);
        val(k) = uu*ww;
    end
    [a ind] = max(val);
    idx = idx(1:ind);
    
    knn = pdist([testdata(:,idx); traindata(:,idx)], 'euclidean');
    knn = squareform(knn);
    knni = knn(1:size(testdata(:,idx),1),size(testdata(:,idx),1)+1:end);

    [y ind] = sort(knni,2);
    nlabels = trainlabels(ind(:,1:numk));
    indices(testind,:) = nlabels;
    distances(testind,:) = y(:,1:numk);
    predlabels(testind) = mode(nlabels,2);
end



weights = 1./distances.^4;
labelvalues = zeros(length(classlabels),length(U));
for i=1:length(classlabels)
    for j=1:size(indices,2)
        labelvalues(i,indices(i,j)) = labelvalues(i,indices(i,j)) + weights(i,j);
    end
end
% labelvalues(:,6) = labelvalues(:,6)*.5;

[a predlabelsW] = max(labelvalues,[],2);


[cc uu dd ww] = conmatrix( classlabels, predlabels);
disp([num2str(round(1000*dd)/10) repmat('   ', [size(cc,1) 1]) num2str(sum(cc,2))]);
disp(num2str([ round(uu*1000)/10 round(ww*1000)/10]));

disp(' ');

[cc uu dd ww] = conmatrix( classlabels, predlabelsW);
disp([num2str(round(1000*dd)/10) repmat('   ', [size(cc,1) 1]) num2str(sum(cc,2))]);
disp(num2str([ round(uu*1000)/10 round(ww*1000)/10]));





save fieldclass.mat










break







if 1==1
% zscore standardization
mu = mean(allfeatures,1);  st = std(allfeatures,[],1);
allfeatures = (allfeatures - repmat( mu, [size(allfeatures,1) 1]))./ repmat( st, [size(allfeatures,1) 1]);

% logistic function, works 10% better than zscoring
% allfeatures = 1./ (1+exp(-allfeatures));

% normalizing
allfeatures = allfeatures./repmat( sqrt(sum(allfeatures(:,:).^2,2)), [1 size(allfeatures,2)]);

idx = 1:size(allfeatures,2);

U = unique(classlabels);
feat = [];
s = zeros(size(U));
for j=1:length(U)
    idx = find(classlabels==U(j));
    feat{j} = allfeatures(idx,:);
    s(j) = length(idx);
end
idx = ml_stepdisc( feat);

idx = idx(1:65);

numk = round(median(s));
knn = pdist(allfeatures(:,idx), 'euclidean');
knn = squareform(knn);

predlabels = zeros(length(classlabels),1);
distances = zeros(length(classlabels),numk);
indices = zeros(length(classlabels),numk);
for i=1:N
    trainind = find(ismember(protlabels,trainidx{i}));
    testind = find(ismember(protlabels,testidx{i}));

    trainlabels = classlabels(trainind);
    
    knni = knn(testind,:);
    knni(:,testind) = [];
    [y ind] = sort(knni,2);
    nlabels = trainlabels(ind(:,1:numk));
    indices(testind,:) = nlabels;
    distances(testind,:) = y(:,1:numk);
    predlabels(testind) = mode(nlabels,2);
end

weights = 1./distances.^8;
labelvalues = zeros(length(classlabels),length(U));
for i=1:length(classlabels)
    for j=1:size(indices,2)
        labelvalues(i,indices(i,j)) = labelvalues(i,indices(i,j)) + weights(i,j);
    end
end
% labelvalues(:,6) = labelvalues(:,6)*.5;

[a predlabels2] = max(labelvalues,[],2);

[cc uu dd ww] = conmatrix( classlabels, predlabels);
disp(num2str([ round(uu*1000)/10 round(ww*1000)/10]));
[cc uu dd ww] = conmatrix( classlabels, predlabels2);
disp(num2str([ round(uu*1000)/10 round(ww*1000)/10]));

e = eye(size(cc));
ccc = cc;
ccc(find(e)) = 0;
w = sum(ccc,1);
w = 1 - (w / sum(w));
w = repmat(w, [ length(classlabels) 1]);

labelvalues2 = labelvalues.*w;

[a predlabels3] = max(labelvalues2,[],2);

[cc uu dd ww] = conmatrix( classlabels, predlabels3);
disp(num2str([ round(uu*1000)/10 round(ww*1000)/10]));

break
disp(classes');
disp(num2str(round(dd*1000)/10));
disp(' ');
disp(num2str(sum(cc,2)'));
disp(' ');
disp(num2str([ round(uu*1000)/10 round(ww*1000)/10]));


break
end

% allfeatures = 1./ (1+exp(-allfeatures));

numk = 1;
predlabels = zeros(length(classlabels),1);
for i=1:N
    trainind = find(ismember(protlabels,trainidx{i}));
    testind = find(ismember(protlabels,testidx{i}));

    traindata = allfeatures(trainind,:);
    testdata = allfeatures(testind,:);
    trainlabels = classlabels(trainind);
    testlabels = classlabels(testind);
    trainprots = protlabels(trainind);
    testprots = protlabels(testind);
    
    
% zscore standardization
mu = mean(traindata,1);  st = std(traindata,[],1);
traindata = (traindata - repmat( mu, [size(traindata,1) 1]))./ repmat( st, [size(traindata,1) 1]);
testdata = (testdata - repmat( mu, [size(testdata,1) 1]))./ repmat( st, [size(testdata,1) 1]);

% normalizing
traindata = traindata./repmat( sqrt(sum(traindata(:,:).^2,2)), [1 size(traindata,2)]);
testdata = testdata./repmat( sqrt(sum(testdata(:,:).^2,2)), [1 size(testdata,2)]);



idx_sda = 1:size(traindata,2);

        knn = pdist([testdata(:,idx_sda); traindata(:,idx_sda)], 'euclidean');
        knn = squareform(knn);
        knn = knn(1:size(testdata(:,idx_sda),1),size(testdata(:,idx_sda),1)+1:end);
        idx = repmat( trainlabels',[size(knn,1) 1]);
        [y ind] = sort(knn,2);
        nlabels = trainlabels(ind(:,1:numk));
        predlabels(testind) = mode(nlabels,2);

    if 1==2
        if 1==1
            U = unique(classlabels);
            feat = [];
            for j=1:length(U)
                feat{j} = traindata(find(trainlabels==U(j)),:);
            end
            idx_sda = ml_stepdisc( feat);
            disp(length(idx_sda));
            disp(length(unique(idx_sda)));
            % idx_sda = unique(idx_sda);
        else
            idx_sda = 1:size(allfeatures,2);
        end
        
        % idx_sda = idx_sda(1:20);
        
        ma = 0;
        for j = 5:10:min(length(idx_sda),95)
            [model_j,options_j,ma_j] = trainClassifier( trainlabels, traindata(:,idx_sda(1:j)), c, g);
            % predict_label = svmpredict( trainlabels, traindata(:,idx_sda(1:j)), model_j, '-b 1');
            % [cc uu dd ma_j] = conmatrix( trainlabels, predict_label);
            if ma_j>ma
                ma = ma_j;
                model = model_j;
                options = options_j;
                opt_j = j;
            end
        end
        opt_j
        [predict_label accuracy weights] = svmpredict( testlabels, testdata(:,idx_sda(1:opt_j)), model, '-b 1');
  
        [a ind] = max(weights,[],2);
        predlabels(testind) = model.Label(ind);
        allweights(testind,model.Label) = weights;
    end
end

[cc uu dd ww] = conmatrix( classlabels, predlabels)

break
classes'
addpath lib/results
[cc uu dd ww] = conmatrix( classlabels, predlabels);

num2str(round(dd*1000)/10)
num2str(sum(cc,2))
num2str([ round(uu*1000)/10 round(ww*1000)/10])

save('fieldclass.mat');
