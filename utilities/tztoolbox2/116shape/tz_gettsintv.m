function its = tz_gettsintv(ts,len)
%TZ_GETTSINTV Calculate intervals for time points. (Unknown)
%   ITS = TZ_GETTSINTV(TS,LEN)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function its = tz_gettsintv(ts,length)
%OVERVIEW
%   Calculate intervals for time points
%PARAMETERS
%   ts - time points
%   len - length of one period
%       if it is empty, the signal is not periodic
%RETURN
%   its - intervals
%DESCRIPTION
%   
%HISTORY
%   27-Mar-2005 Initial write TINGZ
%SEE ALSO
%   

% ints=zeros(size(ts));
its(2:length(ts))=ts(2:end)-ts(1:end-1);
its(1)=ts(1);

if ~isempty(len)
    if len<ts(end)
        error('time points are out of range');
    end
    its(1)=its(1)+(len-ts(end));
end

