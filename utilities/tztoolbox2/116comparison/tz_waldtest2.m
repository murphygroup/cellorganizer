function pvalue = tz_waldtest2(s1,s2)
%TZ_WALDTEST2 Obsolete.
%
%See also ML_WALDTEST2

%function pvalue = tz_waldtest2(s1,s2)
%
%OVERVIEW:
%   wald test
%PARAMETERS:
%   s1 - sample 1
%   s2 - sample 2
%RETURN:
%   pvalue - pvalue
%DESCRIPTION:
%   wald test: (X1bar-X2bar)/var ~ N(0,1)
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   15-DEC-2003 Modified TINGZ

error(tz_genmsg('of','tz_waldtest2','ml_waldtest2'));

mean1=mean(s1);
mean2=mean(s2);
var1=var(s1,1);
var2=var(s2,1);

mean12=mean1-mean2;
var12=var1/size(s1,2)+var2/size(s2,2);

w=mean12/sqrt(var12);

pvalue=2*normcdf(-abs(w),0,1);