function alpha=tz_fitmnom2(y,p1,p2)
%TZ_FITMNOM2 Infer 2-component multinomial mixture.
%   ALPHA = TZ_FITMNOM2(Y,P1,P2)
%   
%   See also

%   19-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function alpha=tz_fitmnom2(y,p1,p2)
%   
%OVERVIEW:
%   Solve multinomial mixture of 2 components
%   model: g(x)=alpha*f1(x)+(1-alpha)*f2(x)
%PARAMETERS:
%   y - data
%   p1 - first component
%   p2 - second component
%RETURN:
%   alpha - coefficient
%DESCRIPTION:
%   maximize the likelihood by Newton's method
%
%HISTORY:
%   15-MAR-2004 Initial write TINGZ
%   04-NOV-2004 Modified TINGZ
%       - add comments

zeroc=find((p1+p2)==0);
y(zeroc)=[];
p1(zeroc)=[];
p2(zeroc)=[];

alpha=0.5;

maxiter=1000;
mine=1e-5;

for i=1:maxiter
    d1=sum(y.*(p1-p2)./(alpha*(p1-p2)+p2));
    d2=-sum(y.*(p1-p2).^2./(alpha*(p1-p2)+p2).^2);
    dalpha=d1/d2;
    
    if alpha-dalpha>1
        alpha=(1+alpha)/2;
    else if alpha-dalpha<0
            alpha=alpha/2;
            
        else
            alpha=alpha-dalpha;
        end
    end
    
    if dalpha<mine
        break;
    end
       
    if i==maxiter
        warning('max iteration reached');
    end
    
end