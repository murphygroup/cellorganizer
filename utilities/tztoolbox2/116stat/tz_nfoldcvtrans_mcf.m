function [avgcm,avgacc]=tz_nfoldcvtrans_mcf(features,n,method,t,transmeth,transel)

%function [avgcm,avgacc]=tz_nfoldcv_mcf(features,n,method,t)
%
%OVERVIEW:
%   n fold cross validation
%PARAMETERS:
%   features - mcf
%   n - n fold
%   method - classification method
%   t - classification parameters
%RETURN:
%   avgcm - average confusion matrix
%   avgacc - average accuracy
%
%HISTORY:
%   ??-???-2004 Initial TINGZ
%   31-OCT-2004 Modified TINGZ
%       - add comments

nclass=length(features);

for i=1:nclass
    norgset(i)=size(features{i},1);
    perm{i}=randperm(norgset(i));
    permfold(i,:)=[0,cumsum(ml_redistrnum(norgset(i),n))];
end

if strcmp(method,'bpnn')
    preprocess=0;
else
    preprocess=1;
end

for k=1:n
    k
    for i=1:nclass
        pf=perm{i};
        testsel=pf(permfold(i,k)+1:permfold(i,k+1));
        testset{i}=features{i}(testsel,:);
        pf(permfold(i,k)+1:permfold(i,k+1))=[];
        trainset{i}=features{i}(pf,:);
    end
    [trainfeats,trainclass,traincellidx]=ml_mcf2combfeats(trainset);
    [testfeats,testclass,testcellidx]=ml_mcf2combfeats(testset);
    [trainfeats,testfeats]=ml_featurenorm(trainfeats,testfeats);
    X=trainfeats(:,transel{1});
%     X=[1./X(:,1),X];
    Y=trainfeats(:,transel{2});
    if ~strcmp(transmeth,'og')
        B=tz_transformfeats(X(traincellidx<=transel{3},:),Y(traincellidx<=transel{3},:),transmeth);
        Yh=tz_evaltransf(testfeats(:,transel{1}),B);
    else
        Yh=testfeats(:,transel{2});
    end
    tg=tz_classify(Yh,Y,trainclass,method,preprocess,t{:});
    [ncvcm{k},pcvcm{k},cvacc(k)]=ml_summaryclassif(testclass,tg);
end

[avgcm,avgacc]=ml_calcavgcm(ncvcm);