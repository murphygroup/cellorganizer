function [pvalue,ts]=tz_kstest2(s1,s2)
%TZ_KSTEST2 Obsolete.
%
%See also ML_KSTEST2

%function pvalue=tz_kstest2(s1,s2)
%OVERVIEW:
%   KS test
%PARAMETERS:
%   s1 - sample 1
%   s2 - sample 2
%RETURN:
%   pvalue - p-value
%   ts - test statistic
%DESCRIPTION
%   just call kstest2 funcion in matlab
%
%HISTORY
%   06-OCT-2003 Initial write TINGZ
%   01-NOV-2004 Modified TINGZ
%       - add ts to return value
error(tz_genmsg('of','tz_kstest2','ml_kstest2');

[h,pvalue,ts]=kstest2(s1,s2);