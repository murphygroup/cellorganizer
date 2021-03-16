function tz_plotgabor(range,sigma,u,phi)
%TZ_PLOTGABOR Plot 2D Gabor function.
%   TZ_PLOTGABOR(RANGE,SIGMA,U,PHI) plots gabor function by specifying
%   several parameters. RANGE is a 2x2 matrix specifying the range of 
%   plot. The first row is the range of X coordinate and the second
%   row is for Y coordinate. SIGMA is a vector of the standard deviations
%   of X and Y, which are SIGMA(1) and SIGMA(2) repectively. U is the
%   requency and PHI is the phase. If PHI is empty, the fourier transform
%   of 0 phage of gabor will be plotted.
%   
%   See also

%   ??-???-2004 Initial write T. Zhao
%   02-NOV-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('Exactly 4 arguments are required')
end

step=(range(:,2)-range(:,1))/100;

[X,Y]=meshgrid(range(1,1):step(1):range(1,2),range(2,1):step(2):range(2,2));

if ~isempty(phi)
    Z=exp(-((X/sigma(1)).^2+(Y/sigma(2)).^2)).*cos(2*pi*u*X+phi);
    
else
    A=2*pi*prod(sigma);
    sigma=sigma/2/pi;
   
    Z=A*(exp(-(((X-u)/sigma(1)).^2+(Y/sigma(2)).^2/2))+ ...
        exp(-(((X+u)/sigma(1)).^2+(Y/sigma(2)).^2/2)));
end

mesh(Z);