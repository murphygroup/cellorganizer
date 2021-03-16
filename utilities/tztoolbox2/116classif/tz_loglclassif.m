function [cm,avgacc]=tz_loglclassif(objcom,objnum,testobj,labels)
%TZ_LOGLCLASSIF Classification by multinomial distribution.
%   CM = TZ_LOGLCLASSIF(OBJCOM,OBJNUM,TESTOBJ,LABELS) returns confusion
%   matrix of classifying data matrix TESTOBJ. OBJCOM is the training
%   data matrix. The rows are cells and columns are numbers of object
%   types. OBJNUM is a vector containing the number of objects in each
%   class. TESTOBJ is the data matrix of testing data, also with cells
%   along rows and object numbers along cloumns. LABELS are true labels
%   of testing data.
%
%   [CM,AVGACC] = TZ_LOGLCLASSIF(OBJCOM,OBJNUM,TESTOBJ,LABELS) also
%   returns average accuracy.

%   14-MAR-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('At least 4 arguments are required')
end

nobj=size(testobj,1);
nclass=size(objcom,1);

logl=tz_testlogl(objcom,objnum,testobj);

[a,b]=max(logl')

cm=zeros(nclass,nclass);

for i=1:nobj
    cm(labels(i),b(i))=cm(labels(i),b(i))+1;
end

avgacc=sum(diag(cm))/sum(cm(:));
cm=cm./repmat(sum(cm,2),1,nclass)*100;