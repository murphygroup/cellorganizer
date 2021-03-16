function [avgcm,avgacc,ncm] = tz_nfoldcv_pslid(features,n,regmeth,t)
%TZ_NFOLDCV_PSLID Obsolete.

%function [avgcm,avgacc,ncvcm] = tz_nfoldcv_pslid(features,n,regmeth,t)
%OVERVIEW
%   
%PARAMETERS
%   features - cell array of feature matrices
%   n - number of cross validation
%   regmeth - classification method
%   t - classification parameters (empty: default parameters)
%RETURN
%   avgcm - average confusion matrix
%   avgacc - average accuracy
%DESCRIPTION
%   
%HISTORY
%   22-Jun-2005 Initial write TINGZ
%SEE ALSO
%   

error(tz_genmsg('of','tz_nfoldcv_pslid','ml_nfoldcv'));

nclass=length(features);

for i=1:nclass
    norgset(i)=size(features{i},1);
    perm{i}=randperm(norgset(i));
    permfold(i,:)=[0,cumsum(tz_redistrnum(norgset(i),n))];
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
    if regmeth==1
        for u=1:max(trainclass)
            if sum(trainclass==u)<=size(trainfeats,2)
                disp('PSLID warning: cross validation error - The number of images must exceed the number features for LDA.');
            end
        end
    end
    regmodel=tz_trainclassif_pslid(trainfeats,trainclass,t,regmeth);
    tg=tz_testclassif_pslid(testfeats,regmodel);
    [ncvcm{k},pcvcm{k},cvacc(k)]=tz_summaryclassif(testclass,tg);
end

[avgcm,avgacc,ncm]=tz_calcavgcm(ncvcm);