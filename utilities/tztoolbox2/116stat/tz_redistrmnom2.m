function [x1,x2,logll]=tz_redistrmnom2(y,p1,p2,N1)

%function [x1,x2,logll]=tz_redistrmnom2(y,p1,p2,N1)
%
%OVERVIEW
%   Fit 2-mixture multinomials when the number for each component is fixed
%PARAMETERS:
%   y - data
%   p1 - first multinomial
%   p2 - second multinomial
%   N1 - number in the first component
%RETURN:
%   x1 - first component
%   x2 - second component
%   logll - log likelihood
%DESCRIPTION
%
%HISTORY:
%   13-MAR-2004 Initial write TINGZ
%   14-MAR-2004 Modified TINGZ
%   05-NOV-2004 Modified TINGZ
%       - add comments

N2=sum(y)-N1;
n=length(p1);
x1=tz_mnornd(N1,p1(1:end-1),1);
while(any(x1>y))
    diff=x1-y;
    [maxdiff,maxpos]=max(diff);
    [mindiff,minpos]=min(diff);
    x1(maxpos)=x1(maxpos)-maxdiff;
    x1(minpos)=x1(minpos)+maxdiff;
end

x2=y-x1;

pdiv=p1./p2;
inc1=x2./(x1+1).*pdiv;
dec1=x1./(x2+1)./pdiv;

[dllinc,posinc]=max(inc1);
tempdec1=dec1;
tempdec1(posinc)=0;
[dlldec,posdec]=max(tempdec1);

while(dllinc*dlldec>1)
    x1(posinc)=x1(posinc)+1;
    x2(posinc)=x2(posinc)-1;
    x1(posdec)=x1(posdec)-1;
    x2(posdec)=x2(posdec)+1;
    inc1(posinc)=x2(posinc)./(x1(posinc)+1).*pdiv(posinc);
    dec1(posdec)=x1(posdec)./(x2(posdec)+1)./pdiv(posdec);

    [dllinc,posinc]=max(inc1);
    tempdec1=dec1;
    tempdec1(posinc)=0;
    [dlldec,posdec]=max(tempdec1);    
end

logll=tz_mnomlogl(x1,p1)+tz_mnomlogl(x2,p2);
% p1(x1==0)=1;
% p2(x2==0)=1;
% 
% for i=1:n
%     logfact1(i)=sum(log(1:x1(i)));
%     logfact2(i)=sum(log(1:x2(i)));
% end
% logll=sum(log(1:N1))+sum(log(1:N2))-sum(logfact1)-sum(logfact2)+sum(x1.*log(p1))+sum(x2.*log(p2));
