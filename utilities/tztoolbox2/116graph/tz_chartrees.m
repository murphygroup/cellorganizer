function [pdeg,treefeat]=tz_chartrees(trees,Vs)
%TZ_CHARTREES Characterize trees.
%   PDEG = TZ_CHARTREES(TREES,VS) returns the vector of percentage of 
%   occurance of node degrees, which are estimated from several trees
%   in the cell array TREES. Vs is a cell array of corresponding node
%   coordinates.
%   
%   [PDEG,TREEFEAT] = TZ_CHARTREES(...) also returns average feature
%   vector of the input trees.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

ntree=length(trees);

degs=[];

for i=1:ntree
    [feats(i,:),name,degree]=tz_calctreefeat(trees{i}',Vs{i});
    degs=[degs,degree];
end

maxdeg=max(degs);
pdeg=zeros(1,maxdeg);

for i=1:length(degs)
    pdeg(degs(i))=pdeg(degs(i))+1;
end

treefeat=mean(feats,1);
pdeg=pdeg/sum(pdeg);


