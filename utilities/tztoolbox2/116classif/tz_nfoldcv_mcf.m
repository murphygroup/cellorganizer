function [avgcm,avgacc]=tz_nfoldcv_mcf(features,n,method,varargin)
%TZ_NFOLDCV_MCF Obsolete.

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
    permfold(i,:)=[0,cumsum(tz_redistrnum(norgset(i),n))];
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
    [trainfeats,trainclass,traincellidx]=tz_mcf2combfeats(trainset);
    [testfeats,testclass,testcellidx]=tz_mcf2combfeats(testset);
    
    tg=tz_classify(testfeats,trainfeats,trainclass,method,preprocess,varargin{:});
    [ncvcm{k},pcvcm{k},cvacc(k)]=tz_summaryclassif(testclass,tg);
end

[avgcm,avgacc]=tz_calcavgcm(ncvcm);