function [pvalue,ts]=tz_chi2test2(s1,s2)
%TZ_CHI2TEST2 Two sample Chi square test.
%   PVALUE = TZ_CHI2TEST2(S1,S2) returns the pvalue of two groups of
%   samples using Chi square test. S1 and S2 are both row vectors and
%   must have the same length. S1 and S2 are supposed to be binning data,
%   so all of their elements should be non-negative integers to get the 
%   meaningful results. The sum of either S1 or S2 should be positive.
%   
%   [PVALUE,TS] = TZ_CHI2TEST2(...) also return the test statistic.
%
%   reference 
%   www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/chi2samp.htm

%   SEP-01-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

if any([s1 s2]<0)
    error('negative number');
end
    
r=sum(s1);
s=sum(s2);

if s==0 | r==0
    error('empty sample');
end

% k1=sqrt(s/r);
% k2=1/k1;

sel=find(s1+s2~=0);
k=length(sel);
s1=s1(sel);
s2=s2(sel);

% ts=sum((k1*s1-k2*s2).^2./(s1+s2));
ts = sum((s*s1-r*s2).^2./(s1+s2)/r/s);  %this is to avoid calc sqrt

if s==r
    c=1;
else
    c=0;
end

if ts==0
    pvalue=1;
else
    pvalue=1-chi2cdf(ts,k-c);
end
