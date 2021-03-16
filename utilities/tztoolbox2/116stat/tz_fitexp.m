function beta=tz_fitexp(x,show,missed)
%TZ_FITEXP Estimate beta in exponential distribution.
%   BETA = TZ_FITEXP(X)
%   
%   BETA = TZ_FITEXP(X,SHOW)
%   
%   BETA = TZ_FITEXP(X,SHOW,MISSED)
%   
%   See also

%   19-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function beta=tz_fitexp(x,show,missed)
%OVERVIEW:
%   Estimate beta in exponential distribution
%PARAMETERS:
%   x - data
%   show - plot the data and fitted curve
%   missed - 
%RETURN:
%   estimated beta
%DESCRIPTION:
%   
%HISTORY:
%   ??-???-???? TINGZ

if(~exist('show','var'))
    show=0;
end

if(~exist('missed','var'))
    missed=0;
end

if missed==0
    beta=mean(x);
else
    beta=mean(x)-1;
end
 
if show~=0
    [n,d] = hist(x,40);
    barplot = bar(d,n/length(x)/(d(2)-d(1)))
    set(barplot,'EdgeColor',[1 1 1]);
    set(barplot,'FaceColor',[0.5 0.5 0.5]);
    
    hold on

    px=tz_ppoints(min(x),max(x),100);
    py=exppdf(px,beta);
    plot(px,py,'r-');
    hold off
end