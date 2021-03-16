function beta = tz_logreg(x,y,n)
%TZ_LOGREG Obsolete. See ML_LOGREG.
%   BETA = TZ_LOGREG(X,Y) returns the maximum likelihood estimates of the
%   parameters BETA in the logistic regression model:
%       P(Y=1|X) = 1/(1+exp(-X*BETA))
%   
%   BETA = TZ_LOGREG(X,Y,N) does the same thing except that the regression
%   model is:
%       Y|X ~ binomial(N,1/(1+exp(-X*BETA))).
%
%   See also

%   27-Apr-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_logreg','ml_logreg'));

if nargin < 2
    error('2 or 3 arguments are required')
end

if nargin<3
    n=1;
end

maxiter=100;
mine=1e-5;

nsample=length(y);
nvar=size(x,2);
beta=zeros(nvar,1);
% sy=y;
% sy(y==1)=1;
% sy(y==0)=1;
sx=1./(1+exp(-x*beta));

% psx1=sx(y==1);
% psx0=sx(y==0);
for i=1:nvar
    delta(i,1)=sum(x(:,i).*(y-n*sx));
    wx=n*x(:,i).*sx.*(1-sx);
    H(i,:)=wx'*x;
end
% llk=sum(log(psx1))+sum(log(psx0));


iH=inv(H);
beta=beta+iH*delta;
llk=0;
for i=1:maxiter
%     beta
%     delta
%     H
    oldbeta=beta;
    sx=1./(1+exp(-x*beta));
    oldllk=llk;
    
    llk=-sum((x(y==0,:)*beta))+sum(log(sx));
    
    if abs(oldllk-llk)<mine | llk==-Inf
        break
    end
    
%     llk=sum(log(psx1))+sum(log(psx0));
    for i=1:nvar
        delta(i,1)=sum(x(:,i).*(y-n*sx));
        wx=n*x(:,i).*sx.*(1-sx);
        H(i,:)=wx'*x;
    end
    
    iH=inv(H);
    beta=beta+iH*delta;
    if mean(abs(beta-oldbeta))<mine
        break
    end
    
end
