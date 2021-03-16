function regmodel = tz_bpnnunmix_clst(training,group,t)
%TZ_BPNNUNMIX_CLST Unknown
%   REGMODEL = TZ_BPNNUNMIX_CLST(TRAINING,GROUP,T)

%   28-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%Unknown
%function y = tz_bpnnunmix_clst(sample,training,group,t)
%OVERVIEW
%   bpnn unmix data by object numbers in each cluster
%PARAMETERS
%   sample - testing data
%   training - training data
%   group - labels of training data
%   t - parameters
%RETURN
%   y - output membership
%DESCRIPTION
%   
%HISTORY
%   29-May-2005 Initial write TINGZ
%SEE ALSO
%   

constidx=tz_constidx(training);
training(:,constidx)=[];
sample(:,constidx)=[];

%fundamental patterns
fp=max(group);

%mixing training size (not counting fundamental patterns)
trainsize=2000;
trainfp=0;
traget = 'wrc';

pmin=1;
rns=[0:3];

%10-fold

for j=1:length(fp)
    trainset=training(group==j,:);
    ntrain(j)=size(trainsel{j},1);
end

%indices of cells for mixing
ids = tz_genmixsel(ntrain,[],trainsize,rns,pmin);

t.prep=1;
t.epochs=300;
t.bp=100;
t.hidden=20;
t.randtrainsel=[];

cellnums=[];

for i=1:size(ids{1},1)
    cellfeats=[];    
    for j=1:length(fp)
        classidx=j;
        if any(cellidx)>0
            cellidx=ids{j}(i,:);
            cellidx=cellidx(cellidx>0);
            cellfeats=[cellfeats;trainset{j}(cellidx,:)];
        else
            cellidx=[];
        end
        
        cellnums(i,j)=length(cellidx);
    end    
    mixcellfeat(i,:)=sum(cellfeats,1);    
end

switch target
case 'frc'
    ytr=tz_normobjcom(cellnums{pfold},0);
case 'num'
    ytr=cellnums{pfold};
case 'wrc'
    ytr=cellnums{pfold}/10;        
end

regmodel = tz_bpnnreg(mixcellfeats,ytr,t);

