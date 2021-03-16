function bias=tz_simmixbias(alpha,B)

%function bias=tz_simmixbias(alpha,B)
%OVERVIEW:
%   generate the random bias for mixture model estimation
%PARAMETERS:
%   alpha - real component weights
%   B - trial times
%RETURN:
%   bias - trial result
%DESCRIPTION:
%
%HISTORY:
%   26-MAY-2004 Initial write TINGZ

k=length(alpha);
w=tz_wrnd(k,B);

for i=1:B
    bias(i)=sum(abs(w(i,:)-alpha));
end