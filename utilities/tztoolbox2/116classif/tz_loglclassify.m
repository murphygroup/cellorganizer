function class = tz_loglclassify(sample,training,group,t)
%TZ_LOGLCLASSIFY Classification by log likelihood.
%   CLASS = TZ_LOGLCLASSIFY(SAMPLE,TRAINING,GROUP,T) classify data by
%   comparing the likelihoods in the classes. The inputs and output are 
%   the same as CLASSIFY except the additional input argument T, which 
%   is a cell array specifying the method of likelihood estimation. 
%   See TZ_ESTLK for more details.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('At least 4 arguments are required')
end

caclass = tz_findclass(group);
nclass = length(caclass);

for i=1:nclass
    tr{i}=training(group==caclass{i}(1),:);
end

for i=1:nclass
    lk(:,i)=tz_estlk(sample,tr{i},1,t{:});
end

nsample=size(sample,1);
for j=1:nsample
    [maxlk,class(j)]=max(lk(j,:));
    class(j)=caclass{class(j)}(1);
end

class=class';

