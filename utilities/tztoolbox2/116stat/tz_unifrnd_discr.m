function sample=tz_unifrnd_discr(nsample,k)
%TZ_UNIFRND_DISCR (Obsolete) Generate discrete uniform numbers.
%   SAMPLE = TZ_OBJFEAT(NSAMPLE,K) returns NSAMPLE integers that are 
%   randomly picked up from the set [1 2 ... K].
%
%   Notice: the function UNIDRND in matlab statistics toolbox does the
%   similiar job. But the input arguments are a bit different. See UNIDRND
%   for more details.
%   
%   See also UNIDRND

%   ??-???-???? Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU   

%% old codes%%
% for i=1:nsample
%     r=randperm(k);
%     sample(i)=r(1);
% end
%%%%%%%%%%%%

sample=floor(rand(1,nsample)*k)+1;

