function [pvalue,ts] = tz_ht2test2(X, Y, pooled)
%TZ_HT2TEST2 Obsolete.
%
%See also ML_HT2TEST2

%[pvalue,ts] = tz_ht2test2(X, Y, pooled)
%OVERVIEW:
%   Hotelling T2 two-sample test
%PARAMETERS:
%   X - sample 1
%   Y - sample 2
%   pooled - pooled or not
%RETURN:
%   pvalue - p-value
%   ts - test statistic
%DESCRIPTION
%   Pooled Hotelling T2 test: 
%       (mean(X)-mean(Y))*inv(Sp*sqrt(1/m+1/n))*(mean(X)-mean(Y))'~(m+n-2)*p/(m+n-p-1)*Fp,n+m-p-1
%   Unpooled Hotelling T2 test: 
%       (mean(X)-mean(Y))*inv(Sx/m+Sy/n)*(mean(X)-mean(Y))'~X2p 
%       (X2 means chi square)
%
%HISTORY
%   ??-???-???? Initial write TINGZ
%   27-MAR-2003 Mofied TINGZ
%   23-NOV-2003 Mofied TINGZ
%   14-DEC-2003 Mofied TINGZ
%   17-FEB-2004 Mofied TINGZ

error(tz_genmsg('of','tz_ht2test2','ml_ht2test2'));
pvalue=1;

if ~exist('pooled','var')
    pooled=1;
end

varX=sum(abs(X(1:end-1,:)-X(2:end,:)));
varY=sum(abs(Y(1:end-1,:)-Y(2:end,:)));
constidx=find(varX+varY==0);

if ~isempty(constidx)
    warning('constant included'); 
    if any(X(1,constidx)~=Y(1,constidx))
        pvalue=0;
        ts=Inf;
    else
        warning(['feature ' num2str(constidx) ' removed'])
        X(:,constidx)=[];
        Y(:,constidx)=[];
    end
end

if pvalue~=0
    if(pooled==0)
        ts=tz_twomaha(X,Y,0,1);
        pvalue=1-chi2cdf(ts,size(X,2));
    else
        [ts,df1,df2]=tz_fvalue(X,Y);
        pvalue=1-fcdf(ts,df1,df2);
    end
end