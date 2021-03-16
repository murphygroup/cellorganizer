function [alpha,loglk,succ]=tz_fitmnomk_fix(y,p,k,weights,minalpha)
%TZ_FITMNOMK_FIX Multinomial mixture decomposition by fixed number of comp.
%   ALPHA = TZ_FITMNOMK_FIX(Y,P,K,WEIGHTS,MINALPHA) returns the mixture
%   coefficients which may generate Y from P.
%   
%   [ALPHA,LOGLK,SUCC] = TZ_FITMNOMK_FIX(...) also returns the log likehood
%   and the flag of fail.
%   
%   See also

%   19-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

ncom=size(p,1);

if isempty(weights)
    weights=ones(1,ncom);
end

[alpha,loglk,succ]=tz_fitmnomk(y,p,minalpha);

if (sum(alpha~=0)<=k)
    return;
end


[salpha,pos]=sort(-alpha.*weights);

pos=pos(1:k);

[subalpha,loglk,succ]=tz_fitmnomk(y,p(pos,:),minalpha);

alpha=zeros(size(alpha));
alpha(pos)=subalpha;
