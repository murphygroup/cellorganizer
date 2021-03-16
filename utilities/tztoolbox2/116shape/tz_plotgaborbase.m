function tz_plotgaborbase(mu,sigma)
%TZ_PLOTGABORBASE Plot Gabor base.
%   TZ_PLOTGABORBASE(MU,SIGMA) plots gabor shape basis with mean MU and
%   deviation SIGMA. The parametric function is:
%        Y(s)=exp(-(s-MU)^2/2/SIGMA^2) * ...
%           [cos(2*pi*(s-MU)/SIGMA);
%            sin(2*pi*(s-MU)/SIGMA)]
%   
%   See also

%   30-Mar-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

s=linspace(mu-.5,mu+.5,10000);
env=exp(-(s-mu).^2/2/sigma^2);
x=env.*cos(2*pi*(s-mu));
y=env.*sin(2*pi*(s-mu));
normx=((-2*pi*sin(2*pi*(s-mu))+(s-mu).*cos(2*pi*(s-mu))/sigma^2).*env).^2;
normy=((2*pi*cos(2*pi*(s-mu))+(s-mu).*sin(2*pi*(s-mu))/sigma^2).*env).^2;
x=x/sum(sqrt(normx+normy)*0.0001);
y=y/sum(sqrt(normx+normy)*0.0001);
% x=x./sqrt(sum(x.^2))*1000;
% y=y./sqrt(sum(y.^2))*1000;
x=(x-min(x));
% y=(y-min(y));
plot(x,y);
axis('square');
