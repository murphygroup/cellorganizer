function [pvalues,ts]=tz_compperm_m(feats1,feats2,labeled,method,n,isplot)
%TZ_COMPPERM_M Under construction.

%TZ_COMPFEATS_M Compare randomly mixed pairs.
%   PVALUES = TZ_COMPFEATS_M(FEATS1,FEATS2,LABELED,METHOD,N,ISPLOT)
%   randomly mix and split FEATS1 and FEATS2 then returns pvalues of
%   these pairs for N trials.
%   
%   [PVALUES,TS] = TZ_COMPFEATS_M(...)
%   
%   See also

%   04-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [pvalues,ts]=tz_compperm_m(feats1,feats2,labeled,method,n,isplot)
%   
%OVERVIEW:
%   do n trials of multivariate testing for random mixed pair
%PARAMETERS:
%   feats1 - feature matrix n1xp
%   feats2 - feature matrix n2xp
%   labeled - labeled or not
%   method - test method
%   n - trial times
%   isplot - qqplot or not
%RETURN:
%   pvalues - pvalues
%   ts - test statistic
%DESCRIPTION:
%   The function can be used to validate the test method
%
%HISTORY:
%   ??-???-???? Initial write TINGZ

error('The function tz_compperm_m is under construction.');

for i=1:n
    [pvalues(i),ts(i)]=tz_compfeats_m(feats1,feats2,labeled,method,'mx');
end

if isplot~=0
    figure
    hh=qqplot(0:(1/(n-1)):1,pvalues);
    delete(hh(2:3));
    hold on
    plot([0,1],[0,1],'r-.')
    axis([0,1,0,1])
    xlabel('uniform(0,1)');
    ylabel('quantile of pvalues');
    hold off
end