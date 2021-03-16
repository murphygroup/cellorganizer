function tz_plothist(x,varargin)
%TZ_PLOTHIST
%   TZ_PLOTHIST(X)
%   
%   See also

%   07-Dec-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

[n,c] = hist(x);

bandwidth = c(2)-c(1);
halfBandwidth = bandwidth/2;

c = c-halfBandwidth;

stairs(c,n,varargin{:});
% ylim([-1 15]);
% plot(xx,yy,varargin);