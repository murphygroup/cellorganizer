function dist=tz_twomaha(X,Y,pooled,zscored)
%TZ_TWOMAHA Obsolete.
%
%See also ML_TWOMAHA

%function dist=tz_twomaha(X,Y,pooled,zscored)
%
%OVERVIEW:
%   calculate mahalanobisis distance 
%PARMATERS:
%   X - sample 1
%   Y - sample 2
%   pooled - pooled or not
%   zscored - zscored or not
%RETURN:
%   dist - mahalanobasis distance
%
%HISTORY:
%   23-NOV-2003 Initial write TINGZ

error(tz_genmsg('of','tz_twomaha','ml_twomaha'));

p=size(X,2);
feats=[X;Y];
m=size(X,1);
n=size(Y,1);

if(zscored~=0)
    feats=zscore(feats);
    X=feats(1:m,:);
    Y=feats((m+1):end,:);
end

if(pooled==0)
    S=cov(X)/m+cov(Y)/n; 
else
    S=(1/m+1/n)*((m-1)*cov(X)+(n-1)*cov(Y))/(m+n-2);
end

mux=mean(X);
muy=mean(Y);

dist=((mux-muy)/S)*(mux-muy)';

