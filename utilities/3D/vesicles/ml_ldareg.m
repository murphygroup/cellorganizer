function regmodel = ml_ldareg(x,y,t)

%ML_LDAREG trains an LDA classifier
%   REGMODEL=ML_LDAREG(X,Y,T) returns a structure restoring 
%   the trained model. T is the structure of
%   parameters. It has the following fields:
%       norm, stop, randtrainsel: see ML_REGRESS
%
%   SEE ALSO ML_EVALLDAREG

% Copyright (C) 2006  Murphy Lab
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

%   29-May-2005 Initial write TINGZ

if (size(y,2)>1)
    error('The target y must have only one column.');
end

ngroup = max(y);
gr=size(y,1);

constidx=[];
for i=1:ngroup
    gx=x(y==i,:);
    constidx=[constidx ml_constidx(gx)];
end

if isempty(constidx)
    prep.featidx=[];
else
    allidx=1:size(x,2);
    allidx(constidx)=[];
    prep.featidx=allidx;
    x(:,constidx)=[];
end

[rx,cx] = size(x);
t.trsize=rx;


if rx ~= gr,
    regmodel.succ=0;
    regmodel.warnmsg='The number of rows in the second and third input arguments must match.';
    error(regmodel.warnmsg); 
end

for k = 1:ngroup
   groupk = x(find(y == k),:);
   [rg,cg]=size(groupk);
   if rg < cg
       regmodel.succ=0;
       regmodel.warnmsg='The number of samples must exceed the number features for LDA.';
       error(regmodel.warnmsg); 
   end
   meanx = mean(groupk);
   [Q,R] = qr(groupk - meanx(ones(rg,1),:),0);
   R = R / sqrt(rg-1);
   trained.rs{k}=R;
   trained.means{k}=meanx;
   trained.prior(k)=rg/rx;
   %tz+ 01-Mar-2007
   s = svd(R);
   trained.logsigma(k) = 2*sum(log(s));
   %tz++
end

regmodel.modelname='lda';
regmodel.type='llk'; %likelihood
regmodel.trained=trained;
regmodel.t=t;
regmodel.prep=prep;
regmodel.postp.ctg=1
