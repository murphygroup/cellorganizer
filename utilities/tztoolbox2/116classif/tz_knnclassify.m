function class = tz_knnclassify(sample,training,group,k)
%TZ_KNNCLASSIFY KNN classification.
%   CLASS = TZ_KNNCLASSIFY(SAMPLE,TRAINING,GROUP,K) is similiar with
%   CLASSIFY. But it has one more argument K to determine the number
%   of nearest neighbors.
%   
%   See also CLASSIFY

%   28-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

class = tz_knn(training, sample, tz_label2post(group), k);