function pvalue=tz_chisqtest_multinom(x,p)
%TZ_CHISQTEST_MULTINOM One sample Chi square test.
%   PVALUE = TZ_CHISQTEST_MULTINOM(X,P)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function pvalue=tz_chisqtest_multinom(x,p)
%   
%OVERVIEW:
%   chi square test
%PARAMETERS:
%   x - sample
%   p - multinomial parameters
%RETURN:
%   pvalue - p-value
%DESCRIPTIOIN
%
%HISTORY:
%   10-MAR-2003 nitial write TINGZ

N=sum(x);
ex=p*N;

if(any(ex==0))
    impidx=find(ex==0);
    if(any(x(impidx)~=0))
        pvalue=0;
        return;
    end
    
    
    x(impidx)=[];
    ex(impidx)=[];
end
deg=length(x)-1;

chisq=sum((x-ex).^2./ex);
pvalue=1-chi2cdf(chisq,deg);