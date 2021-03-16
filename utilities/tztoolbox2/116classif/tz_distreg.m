function regmodel = tz_distreg(x,y,t)
%TZ_DISTREG Train a nearest center model.
%   REGMODEL = TZ_DISTREG(X,Y) returns a model that is trained on X. The
%   model is to desribe the center of each class. It is only for categorized
%   data.
%   
%   REGMODEL = TZ_DISTREG(X,Y,T) also let the user specify parameters.
%   There is no additional parameters currently except t.norm. In the 
%   future version different distance functions will be supported.
%   
%   See also ML_REGRESS

%   21-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('2 or 3 arguments are required')
end

if (size(y,2)>1)
    error('The target y must have only one column.');
end

if ~exist('t','var')
    t.norm=1;
    t.distfun='euc';
end

if t.norm==1
    [x,zmean,zsdev]=ml_zscore(x);
    prep.zscore=1;
    prep.zmean=zmean;
    prep.zsdev=zsdev;
end

nclass = max(y);
for k=1:nclass
    centers(k,:) = sum(x(y==k,:),1);
end

regmodel=struct('modelname','NCC','modeltype',...
    'DIST','trained',centers,'t',t,'prep',prep);
    
regmodel.postp.ctg=1;
   