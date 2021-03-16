function constidx=tz_constidx(X)
%TZ_CONSTIDX Obsolete.
%
%See also ML_CONSTIDX

%function idx=tz_constidx(feats)

error(tz_genmsg('of','tz_constidx','ml_constidx'));

varX=sum(abs(X(1:end-1,:)-X(2:end,:)),1);

constidx=find(varX==0);