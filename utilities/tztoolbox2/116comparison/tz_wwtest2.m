function [pvalue,ts]=tz_wwtest2(X1,X2,s)
%TZ_WWTEST2 Obsolete.
%
%See also ML_WWTEST2

%function [pvalue,ts]=tz_wwtest2(X1,X2,s)
%OVERVIEW:
%   Multivariate Wald-Wolfowitz test
%PARAMETERS:
%   X1 - sample 1
%   X2 - sample 2
%   s - distance functioin
%RETURN:
%   pvalue - p-value
%   ts - test statistic
%DESCRIPTION:
%   FR test. Jerome H. Friedman and Lawrence C. Rafsky, 1979
%HISTORY:
%   05-NOV-2003 Initial write TINGZ
%   17-FEB-2004 Modified TINGZ
%   30-JUL-2004 Modifien

error(tz_genmsg('of','tz_wwtest2','ml_wwtest2'));

% if ~exist('s','var')
%     s='eu';
% end
%     
% X=zscore([X1;X2]);
% 
% m=size(X1,1);
% n=size(X2,1);
% N=m+n;
% v=tz_mstree(X,s,m);
% 
% R=sum((v(1,:)<=m) & (v(2,:)>m))+sum((v(1,:)>m) & (v(2,:)<=m))+1;
% Er=2*m*n/N+1;
% C=0;
% 
% for i=1:N
%     deg=sum(v(1,:)==i)+sum(v(2,:)==i);
%     C=C+deg*(deg-1)/2;
% end
% 
% Vr=2*m*n*((2*m*n-N)/N+(C-N+2)*(N*(N-1)-4*m*n+2)/(N-2)/(N-3))/N/(N-1);
% W=(R-Er)/sqrt(Vr);
% ts=R;
% pvalue=normcdf(W,0,1);
msg = tz_genmsg('of','tz_wwtest2','ml_wwtest2');
error(msg);
