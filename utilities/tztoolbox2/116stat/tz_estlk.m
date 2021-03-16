function lk=tz_estlk(Y,X,method,islog,t)
%TZ_ESTLK Likelihood estimation.
%   LK = TZ_ESTLK(Y,X,METHOD,0,T) returns likelihood of Y based on the
%   density estimated from X. See TZ_TRAINLK for details about METHOD and
%   T.
%
%   LK = TZ_ESTLK(Y,X,METHOD,1,T) returns log-likehood.
%     
%   See also TZ_TRAINLK, TZ_EVALK

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 5
    error('Exactly 5 arguments are required')
end

params=tz_trainlk(X,method,t);
lk=tz_evalk(Y,params,islog);