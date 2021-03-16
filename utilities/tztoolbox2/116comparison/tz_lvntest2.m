function [pvalue,ts]=tz_lvntest2(X,Y)
%TZ_LVNTEST2 Levene's test for variances.
%   PVALUE = TZ_LVNTEST2(X1,X2) returns the p-value of the Levene's test
%   on X1 and X2. This is to test if the variances of X1 and X2 are the
%   same.
%
%   [PVALUE,TS] = TZ_LVNTEST2(...) also returns the test statistic.

%   23-NOV-2003 Initial write TT. Zhao
%   01-NOV-2004 Modified T. Zhao
%       - add ts in return
%       - fix calling old version function
%   Copyright (c) Murphy Lab, Carnegie Mellon University

X=zscore(X);
m=size(X,1);
Y=zscore(Y);
n=size(Y,1);

mux=mean(X);
muy=mean(Y);

X=abs(X-repmat(mux,m,1));
Y=abs(Y-repmat(muy,n,1));

[pvalue,ts]=ml_ht2test2(X,Y);