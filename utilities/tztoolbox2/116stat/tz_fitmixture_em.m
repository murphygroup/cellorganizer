function alphas=tz_fitmixture_em(y,pdfs,ts)
%TZ_FITMIXTURE_EM Fit mixture model by EM algorithm
%   ALPHAS = TZ_FITMIXTURE_EM(Y,PDFS,TS)
%   
%   See also

%   19-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function alphas=tz_fitmixture_em(y,pdfs,ts)
%
%OVERVIEW:
%   fit mixture model by em algorithm
%PARAMETERS:
%   y - data
%   pdfs - pdf functions
%   ts - parameters of the pdfs
%RETURN:
%   alphas - coeffiecients
%DESCRIPTION:
%   MODEL: f(y)=a1*f1(x)+...+ak*fk(x); a1+...+ak=1
%   EM iteration:
%       aj(m+1)=sum(P(y|zj=1)/P(y))/n
%       see trainmixture.doc by tingz
%HISTORY:
%   16-NOV-2004 Initial write TINGZ

k=length(ts);
n=length(y);

%initialize alphas
alphas(1,:)=ones(1,k)/k;

%calculate probability matrix
for i=1:k
    if length(pdfs)==1)
        pdfi=pdfs{1};
    else
        pdfi=pdfs{i};
    end
    fy(i,:)=eval([pdfi '(' '''y,''' ts{i} ')'];
end

niter=1000;
minerr=0.0001;

for i=1:k-1
    alphas(2,i)=alphas(1,i)*sum((fy(i,:)./alphas(1,:)*fy)/n;
end

while max(abs(alphas(2,:)-alphas(1,:)))>minerr
    alphas(1,:)=alphas(2,:);
    for i=1:k-1
        alphas(2,i)=alphas(1,i)*sum((fy(i,:)./alphas(1,:)*fy)/n;
    end
end



