function pvalue=tz_chisqtest2_multinom(x1,x2)

%UUU
%function pvalue=tz_chisqtest2_multinom(x1,x2)
%OVERVIEW:
%   two-sample chisq test
%PARAMETERS:
%   x1 - sample 1
%   x2 - sample 2
%RETURN:
%   pvalue -  p-value
%
%HISTORY:
%   10-MAR-2003 Initial write TINGZ

x=x1;
N=sum(x);
%p=(tz_normobjcom(x1)+tz_normobjcom(x2))/2;
x3=x2(find(x2~=0));
minx=min(x3);

x4=x2;
x4=x4*2;
x4(find(x4==0))=1;
x=x+(x4==0)*0.5;

p=tz_normobjcom(x4);
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