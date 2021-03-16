function alpha=tz_solvemix_fix(y,p,k,weights)
%TZ_SOLVEMIX_FIX Unknown.
%   ALPHA = TZ_SOLVEMIX_FIX(Y,P,K,WEIGHTS)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [alpha,loglk,succ]=tz_fitmnomk_fix(y,p,k,weights)

ncom=size(p,1);

if isempty(weights)
    weights=ones(1,ncom);
end

alpha=tz_solvemix(y,p);

if (sum(alpha~=0)<=k)
    return;
end


[salpha,pos]=sort(-alpha.*weights);

pos=pos(1:k);

subalpha=tz_solvemix(y,p(pos,:));

alpha=zeros(size(alpha));
alpha(pos)=subalpha;
