function n=tz_redistrnum(N,k)
%TZ_REDISTRNUM Obsolete.
%
%See also ML_REDISTRNUM

%function n=tz_redistrnum(N,k)
%
%OVERVIEW:
%   evenly separate a number into several numbers
%PARAMETERS:
%   N - a number
%   k - how many small numbers
%RETURN:
%   n - a vector
%DESCRIPTION:
%   e.g. N=10,k=4-->[3,2,2,3]
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   05-NOV-2004 Modified TINGZ
%       - add comments
%   22-MAR-2005 Modified TINGZ
%       - debug k==1
%   

error(tz_genmsg('of','tz_redistrnum','ml_redistrnum'));

if k==1
    n=N;
    return
end

avgn=floor(N/k);
for i=1:k
    n(i)=avgn;
end

remain=N-sum(n);
k=1;
while remain>0
    n(k)=n(k)+1;
    remain=remain-1;
    k=k+1;
end

% 
% n(k)=N-sum(n(1:k-1));
% 
% [maxn,maxpos]=max(n);
% [minn,minpos]=min(n);
% 
% while(maxn-minn>1)
%     n(maxpos)=n(maxpos)-1;
%     n(minpos)=n(minpos)+1;
%     [maxn,maxpos]=max(n);
%     [minn,minpos]=min(n);
% end