function params = tz_traindiriclk(x)
%TZ_TRAINDIRICLK Train a dirichlet distribution
%   PARAMS = TZ_TRAINDIRICLK(X)
%   
%   See also

%   21-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end