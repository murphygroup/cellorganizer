function class = tz_nnobjclassify(sample,training,group,weighted)
%TZ_NNOBJCLASSIFY Classification based on object neighbors.
%   CLASS = TZ_NNOBJCLASSIFY(SAMPLE,TRAINING,GROUP,WEIGHTED) classifies
%   data by comparing objects between cells. It searches the nearest 
%   neighbor of each object in a cell, which is an element of the cell
%   array SAMPLE, and counts the class labels of the nearest neighbors. 
%   The cell is classified into a class with most objects being the
%   nearest neighbors of the testing cell objects. TRAINING is a matrix
%   containg all traing object features. WEIGHTED sepcifies the 
%   weighting method:
%       1 - use the last column in elements of SAMPLE
%       0 - equal weights
%       -1 - log likelihood       
%   CLASS and GROUP are vectors and have the same meaning as those in
%   CLASSIFY.

%   AUG-26-2004 Initial write T. Zhao
%   SEP-06-2004 Modified T. Zhao
%       - under construction
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('3 or 4 arguments are required')
end

if ~exist('weighted','var')
    weighted=0;
end

nclass=max(group);
ntrainobj=size(training,1);

%check if training set has cell labels
if size(training,2)-size(sample{1},2)==1
    traincellidx=training(:,end);
    training(:,end)=[];
    
    %build histogram for obj numbers of each class
    for i=1:nclass
        ncell=max(traincellidx(group==i));
        for j=1:ncell
            nums(j,1)=sum(group==i & traincellidx==j);
        end
        trainnums{i}=[nums,ones(size(nums,1),1)];
        alltrainnums(i)=sum(nums);
    end
end

nnp=nn_prepare(training);

for i=1:length(sample)
    switch weigthed
    case 1
        testing=sample{i}(:,1:end-1);
        weights=sample{i}(:,end);
    case 0
        testing=sample{i};
        weights=ones(size(sample{i},1),1);
    case -1
        testing=sample{i};
    end
    
    [index, distance] = nn_search(training, nnp, testing, 1);
    
    post=zeros(1,max(group));
    
    for j=1:length(index)
        if weighted==-1
            weight=1/alltrainnums(group(index(j)));
        else
            weight=weights(j);
        end
        post(group(index(j)))=post(group(index(j)))+weights(j); %+weight?
    end
    
    %log density for the objnum of testing set
    if weigthed==-1
        nobj=size(sample{i},1);
        for j=1:nclass
            lognumden(j)=log(tz_discrdens(trainnums{j},50,nobj));
            
        end
        post=post+lognumden;
    end
    
    [ignore,class(i)]=max(post);
end

class=class';