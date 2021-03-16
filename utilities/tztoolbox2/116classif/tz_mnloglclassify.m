function class = tz_mnloglclassify(sample,training,group,t)
%TZ_MLLOGLCLASSIFY Classification based on multinomial models. 
%   CLASS = TZ_MLLOGLCLASSIFY(SAMPLE,TRAINING,GROUP,T) classifies data
%   based on multinomial distribution. TRAINING is the matrix containg
%   object numbers. TRAINING(i,j) is the number objects in cell i of
%   type j. SAMPLE is the testing data and has the same meaning as that
%   of TRAINING. T is an integer specifying how to get object numbers
%   in each class:
%       >0 - the Tth column from the TRAINING
%       0 - sum the rows of TRAINING 
%       -1 - no object numbers
%   GROUP and CLASS are vectors of class labels. See CLASSIFY for more
%   details.

%   ???-??-???? Initial write T. Zhao
%   APR-18-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin < 4
    error('At least 4 arguments are required')
end

if t>0
    combobjnum=training(:,t);
    training(:,t)=[];
    testnum=sample(:,t);
    sample(:,t)=[];
else
    testnum=0;
end

trainset=ml_findclass([training,group]);

for i=1:length(trainset)
    clabel=trainset{i}(1,end-1);
    trainset{i}=trainset{i}(:,1:(end-2));
    if t>0
        objnum{i}=[combobjnum(group==clabel),ones(sum(group==clabel),1)];
    end
end
    
[objsq,objnum2]=tz_trainobjcom(trainset,[]);

if t==0
    objnum=objnum2;
end

if t==-1
    objnum=[];
end

logl=tz_testlogl(objsq,objnum,sample,{},testnum);

[a,class]=max(logl');
class=class';
