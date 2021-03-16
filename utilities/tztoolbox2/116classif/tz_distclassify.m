function class=tz_distclassify(sample,training,group)
%TZ_DISTCLASSIFY Nearest centroid classifier
%   CLASS = TZ_DISTCLASSIFY(SAMPLE,TRAINING,GROUP) does the same thing
%   as CLASSIFY but using different classification method, which is
%   the nearest centroid classifier (NCC) here.
%   
%   See also CLASSIFY

%   ??-???-???? Initial write T. Zhao
%   19-APR-2004 Modified T. Zhao
%       - standarize
%   16-MAY-2004 Modified T. Zhao
%        - correct the centers order
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('At least 3 arguments are required')
end

caclass=ml_findclass(group);
nclass=length(caclass);

for i=1:nclass
    centers(i,:)=mean(training(group==caclass{i}(1),:),1);
end

nsample=size(sample,1);

for i=1:nsample
    for j=1:nclass
        dist(j)=sqrt(sum((sample(i,:)-centers(j,:)).^2));
    end
    [m,class(i)]=min(dist);
end

class=class';