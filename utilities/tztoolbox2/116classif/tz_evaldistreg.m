function y = tz_evaldistreg(x,regmodel)
%TZ_EVALDISTREG
%   Y = TZ_EVALDISTREG(X,REGMODEL)
%   
%   See also

%   21-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

if isfield(regmodel,'prep')
    if isfield(regmodel.prep,'featidx')
        if ~isempty(regmodel.prep.featidx)
            x=x(:,regmodel.prep.featidx);
        end
    end
    if regmodel.prep.zscore
        x=ml_zscore(x,regmodel.prep.zmean,regmodel.prep.zsdev);
    end
end

switch regmodel.t.distfun
case 'euc'
    if ~isfield(regmodel.t,'weights')
        d = dist2(x,regmodel.trained);
        [mindist,y] = min(d,[],2);
    end
otherwise
    error('Invalid distance function in the model.');
end
