function y = tz_evalreg(x,model)
%TZ_EVALREG Evaluate regression model.
%   Y = TZ_EVALREG(X,MODEL)
%   
%   See also

%   10-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end