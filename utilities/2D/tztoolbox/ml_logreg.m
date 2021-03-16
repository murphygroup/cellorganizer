function beta = ml_logreg(x,y,n)
%ML_LOGREG Logistic regression.
%   BETA = ML_LOGREG(X,Y) returns the maximum likelihood estimates of the
%   parameters BETA in the logistic regression model:
%       P(Y=1|X) = 1/(1+exp(-X*BETA))
%   
%   BETA = ML_LOGREG(X,Y,N) does the same thing except that the regression
%   model is:
%       Y|X ~ binomial(N,1/(1+exp(-X*BETA))).
%
%   See also

%   27-Apr-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% Copyright (C) 2007  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu
%
% 12/10/2020 R.F.Murphy use \ operator instead of matrix inversion
%                       (faster and reduces singular matrix errors)

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

% replace inv with \
%iH=inv(H);
%beta=beta+iH*delta;
beta = beta + (H \ delta);
llk=0;

%for testing
% llks = [];
%%

for i=1:maxiter
%     beta
%     delta
%     H
      
    oldbeta=beta;
    sx=1./(1+exp(-x*beta));
    oldllk=llk;
    
    llk=-sum((x(y==0,:)*beta))+sum(log(sx));
    
%for testing
%     llks(i) = llk;
%     plot(llks)
%     drawnow
%%

    if abs(oldllk-llk)<mine | llk==-Inf
        break
    end
    
%     llk=sum(log(psx1))+sum(log(psx0));
    for i=1:nvar
        delta(i,1)=sum(x(:,i).*(y-n*sx));
        wx=n*x(:,i).*sx.*(1-sx);
        H(i,:)=wx'*x;
    end
    
% replace inv with \
    %iH=inv(H);
    %beta=beta+iH*delta;
    beta = beta + (H \ delta);
    if mean(abs(beta-oldbeta))<mine
        break
    end
    
end
