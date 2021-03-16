function gsfilt = tz_gsfilter(sigma,range)
%TZ_GSFILTER Create gaussian filter.
%   IMG2 = TZ_GSFILTER(sigma)
%   
%   See also

%   04-Oct-2005 Initial write T. Zhao

if nargin < 2
    error('At least 2 arguments are required')
end

if nargin < 2
    range = 3*sigma;
end
range = round(range);

gsfilt = normpdf(-range:range,0,sigma);

