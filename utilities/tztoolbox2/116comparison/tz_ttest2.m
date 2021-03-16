function [pvalue,ts] = tz_ttest2(s1,s2)
%TZ_TTEST2 Obsolete.
%
%See also ML_TTEST2

%function [pvalue,ts] = tz_ttest2(s1,s2)
%OVERVIEW:
%   two-sample ttest
%PARAMETERS:
%   s1 - sample 1
%   s2 - sample 2
%RETURN:
%   pvalue - p-value
%DESCRIPTION
%   t test, call ttest2 function but only return pvalue and ts
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   05-NOV-2004 Modified TINGZ
%       - add return value ts

error(tz_genmsg('of','tz_ttest2','ml_ttest2'));

[a,pvalue,c,ts]=ttest2(s1,s2);
    
