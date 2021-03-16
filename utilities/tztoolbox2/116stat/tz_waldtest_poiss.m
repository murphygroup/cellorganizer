function pvalue = tz_waldtest_poiss(s,lamda)

%function pvalue = tz_waldtest_poiss(s,lamda)
%OVERVIEW:
%   wald test for poisson distribution
%PARAMETERS:
%   s - sample
%   lamda - H0
%RETURN:
%   pvalue - p-value
%
%DESCRIPTION:
%   
%HISTORY
%   ??-???-???? Initial write TINGZ
%   15-DEC-2003 Modified TINGZ

act=lamda;
est=mean(s);
se=sqrt(est/size(s,2));

w=(est-act)/se;

pvalue=2*normcdf(-abs(w),0,1);