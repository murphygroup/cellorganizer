function logl=tz_mnomlogl(x,p,pN)

%logl=tz_mnomlogl(x,p,pN)
%OVERVIEW:
%   Calculate log likelihood for multinomial data
%PARAMETERS:
%   x - data
%   p - pramaters of the model
%   pN - pmf of object numbers
%RETURN:
%   logl - log likelihood
%DESCRIPTION:
%
%HISTORY:
%   15-MAR-2004 Initial write TINGZ
%   30-MAR-2004 Modified TINGZ
%       add calculation likelihood for mixture multinomials
%   05-NOV-2004 Modified TINGZ
%       add comments

if ~exist('pN','var')
    pN=ones(size(x,1),1);
end

y=x;
py=p;
nmix=size(x,1);
logl=0;
for k=1:nmix
    x=y(k,:);
    p=py(k,:);
    n=length(x);
    N=sum(x);
   
    if N==0
        continue;
    end
    
    if pN(k)==0
        logl=-Inf;
        break;
    else
        lpN=log(pN(k));
        
        p(x==0)=1;
        
        for i=1:n
            logfact(i)=sum(log(1:x(i)));
        end
        
        logl=logl+sum(log(1:N))-sum(logfact)+sum(x.*log(p))+lpN;
    end
end