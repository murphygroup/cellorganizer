zfunction sp = tz_feat2sp(knots,coefs)
%TZ_FEAT2SP Obsolete See ML_FEAT2SP.
%   SP = TZ_FEAT2SP(KNOTS,COEFS) returns a B-spline structure with internal
%   nodes KNOTS and coefficients COEFS.
%   
%   See also

%   28-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_feat2sp','ml_feat2sp'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

sp.form = 'B-';
sp.coefs = coefs;
sp.number = length(coefs);
sp.order = length(coefs)-length(knots);
sp.dim = 1;
sp.knots = [zeros(1,sp.order) knots ones(1,sp.order)];
