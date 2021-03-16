function sigma=tz_wmoment(data,order)

%function sigma=tz_wmoment(data,order)
%
%OVERVIEW:
%   weighted central moments of all orders.
%PARAMETERS:
%   data - nx2 [value,weight]
%   order - order of the moment
%RETURN:
%   sigma - moment
%DESCRIPTION:
%
%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   03-NOV-2004 Modified TINGZ
%       - add comments

error(tz_genmsg('of','tz_wmoment','ml_wmoment'));

weights=data(:,2);
data=data(:,1);


[m,n] = size(order);
if max(m,n) ~= 1
    error('Requires a scalar second argument.');
end
if  (order - floor(order)) > 0 | order < 1
    error('Requires a positive integer second argument.');
end

sw=sum(weights);
mu = sum(data.*weights)/sw;
   
if order == 1
   sigma = mu;
else
   sigma = sum(((data - mu).^order).*weights)/sw;
end  
