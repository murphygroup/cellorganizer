function g = tz_gaborbank(nscale,norient,freqrange,index,varargin)
%TZ_GABORBANK Babor filter bank.
%   G = TZ_GABORBANK(NSCALE,NORIENT,[UL UH],[M N],[XSIZE YSIZE]) returns
%   a gobor filter from a gabor bank with the index M,N. 
%   
%   See also TZ_GABOR

%   03-Oct-2005 Initial write T. Zhao

if nargin < 4
    error('Exactly 4 arguments are required')
end

Ul = freqrange(1);
Uh = freqrange(2);
S = nscale;
K = norient;
m = index(1);
n = index(2);
W = 0.52;

a = (Uh/Ul)^(1/(S-1));
alpha = (a+1)/(a-1);
sigmax = (1/2/pi)*alpha*( sqrt( 2*log(2)/Uh ) );
sigmay = 1/( tan(pi/2/K)*sqrt( (alpha^2-1)/(sigmax^2) ) );
theta = n*pi/K;

g = tz_gabor([sigmax,sigmay],W,a^m,theta,varargin{:});
