function caclass=ml_findclass(class,ispint)
%ML_FINDCLASS Find classes in a vector or matrix.
%   CACLASS = ML_FINDCLASS(CLASS) separates CLASS into a cell array
%   class by class and returns the cell array. CLASS should be a vector of
%   class labels or a matrix with class labels in the last column. In each
%   matrix of the returned cell array, the last column contains sample
%   labels and the 2nd last column contains class labels. The class labels
%   must be positive integers.
%   
%   CACLASS = ML_FINDCLASS(CLASS,ISPINT) also lets the user specify the
%   type of class label. If the class labels are positive integers,
%   we can set ISPINT to 1 so that the function can be implemented in a
%   faster way. If not, it should be 0, or the function will not be
%   implemented properly.

%   APR-18-2004  Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('ispint','var')
    ispint=1;
end

if ispint==0
  
    caclass={};
    caclass{1}=[class(1,:),1];
    nsample=size(class,1);
    assigned=0;
    
    for i=2:nsample
        assigned=0;
        for j=1:length(caclass)
            if (class(i,end)==caclass{j}(1,end-1))
                caclass{j}=[caclass{j};[class(i,:),i]];
                assigned=1;
            end
        end
        
        if assigned==0
            caclass{length(caclass)+1}=[class(i,:),i];
        end
    end
else
    label=class(:,end);
    post=ml_label2post(label);
    nclass=max(label);
    if any(sum(post,1)==0)
        post=ml_stdclst(post);
        label=ml_post2label(post);
        nclass=max(label);
    end
    samples=1:size(class,1);
    
    %class label and sample indices
    for i=1:nclass
        caclass{i}=[class(label==i,:),samples(label==i)'];
    end
end