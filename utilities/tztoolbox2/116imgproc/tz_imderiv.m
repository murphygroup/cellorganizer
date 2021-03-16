function dimg = tz_imderiv(img,direction)
%TZ_IMDERIV Calculate partial derivative of an image.
%   DIMG = TZ_IMDERIV(IMG,DIRECTION) returns a partial derivative of the
%   image IMG. DIRECTION is used to specify the derivative direction:
%       'x' or 'X' - x direction
%       'y' or 'Y' - y direction
%   The returned value DIMG is a double matrix.
%
%   See also

%   24-Oct-2005 Initial write T. Zhao

if nargin < 2
    error('Exactly 2 arguments are required')
end

img = double(img);

switch direction
    case {'x','X'}
        derivFilter = [-1; 1; 0];
    case {'y','Y'}
        derivFilter = [-1 1 0];
    otherwise
        error('Invalid direction');
end

dimg = imfilter(img,derivFilter,'same');