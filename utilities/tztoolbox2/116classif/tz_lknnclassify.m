function class=tz_lknnclassify(sample,training,group,k)
%TZ_KNNCLASSIFY KNN classification for large samples.
%   CLASS = TZ_KNNCLASSIFY(SAMPLE,TRAINING,GROUP,K) does the same job
%   as TZ_KNNCLASSIFY. But it is designed for large size of samples.
%   
%   See also TZ_KNNCLASSIFY

%   ??-???-???? Initial write T. Zhao
%   04-NOV-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University


% nsample=size(sample,1);
% 
% for i=1:nsample
%     class(i)=tz_knnclassify(sample(i,:),training,group,k);
%     if mod(i,100)==0
%         display([num2str(i) '/' num2str(nsample) ' ']);
%     end
% end
%     
% class=class';          

[index, distance] = nn_search(training, nn_prepare(training), sample, k);

nsample=size(sample,1);

for i=1:nsample
    t=tz_label2post(group(index(i,:)));
    sumt=sum(t,1);
    [tmp,class(i)]=max(sumt+0.5*rand(size(sumt)));
end
class=class';