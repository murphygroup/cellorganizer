function pvalue = tz_waldtest2_binom(x1,n1,x2,n2)
%TZ_WALDTEST2_BINOM Wald test for binomial distribution.
%   PVALUE = TZ_WALDTEST2_BINOM(X1,N1,X2,N2) returns the p-value of the
%   testing on two binomial data X1 and X2 if p1=p2, where
%   X1 ~ binomial(N1,p1), X2 ~ binomial(N2,p2).
%   
%   See also ML_WALDTEST2

%   ??-???-????  Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

p1=x1/n1;
p2=x2/n2;
se=sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2);

delta=p1-p2;

w=delta/se;

pvalue=2*normcdf(-abs(w),0,1);