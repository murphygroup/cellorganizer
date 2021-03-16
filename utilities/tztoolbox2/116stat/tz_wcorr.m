function c = tz_wcorr(x,y,ws)
%TZ_WCORR Weigthed correlation.
%   C = TZ_WCORR(X,Y,WS) returns the weighted correlation between two
%   vectors X and Y. WS is a vector of weights. X, Y and WS must have the
%   same size.
%   
%   See also

%   18-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 3
    error('Exactly 3 arguments are required')
end

mux = ml_wmoment(x,ws,1);
muy = ml_wmoment(y,ws,1);

c = sum((x-mux).*(y-muy).*ws)/sum(ws);
