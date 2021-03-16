function g = tz_gabor(sigmas,bandwidth,scale,theta,range)
%TZ_GABOR Construct a gabor filter
%   G = TZ_GABOR([SIGMAX SIGMAY],BANDWIDTH) returns a 2D gabor
%   filter with the variance SIGMAX at X direction and SIGMAY at Y
%   direction and the requency bandwidth BANDWIDTH.
%   
%   G = TZ_GABOR([SIGMAX SIGMAY],BANDWIDTH,SCALE,THETA,RANGE) also
%   specifies scale, orientation and range.
%
%   See also

%   03-Oct-2005 Initial write T. Zhao

if nargin < 2
    error('At least 2 argument is required')
end

if nargin < 3
    scale = 1;
end

if nargin <5
    range = round(3*sigmas*scale);
    range = [max(range) max(range)];
end

coords = tz_imcoords(range*2+1,scale,(-range-1)/scale);

if nargin == 4
    coords = tz_rotate(coords,theta);
end

g = (1/pi/prod(sigmas))*exp(-(1/2)*(coords(1,:).^2/sigmas(1)^2 ...
    +coords(2,:).^2/sigmas(2)^2)+2*i*pi*bandwidth*coords(1,:))/scale;

g = reshape(g,range(1)*2+1,range(2)*2+1);