function m = tz_wmoment(x,ws,order)
%TZ_WMOMENT Obsolete.
%
%See also ML_WMOMENT

%function m = tz_wmoment(x,ws,order)
%OVERVIEW
%   calculate weighted central moments
%PARAMETERS
%   x - data
%   ws - weights
%   order - 
%RETURN
%   m - moment
%DESCRIPTION
%   the first order will return moment not central moment
%HISTORY
%   24-Apr-2005 Initial write TINGZ
%SEE ALSO
%   

error(tz_genmsg('of','tz_wmoment','ml_wmoment'));

if order==1
    m=sum(x.*ws)/sum(ws);
    return
end

mu=tz_wmoment(x,ws,1);
m=sum((x-mu).^order.*ws)/sum(ws);