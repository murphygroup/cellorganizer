function [alpha,besta,sola]=tz_fitmnom2_ga(y,p1,p2)
%TZ_FITMNOM2_GA Infer 2-component mixture model by gradient ascent.
%   ALPHA = TZ_FITMNOM2_GA(Y,P1,P2)
%   
%   [ALPHA,BESTA,SOLA] = TZ_FITMNOM2_GA(...)
%   
%   See also

%   19-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [alpha,besta,sola]=tz_fitmnom2_ga(y,p1,p2)
%
%OVERVIEW:
%   decompose 2-mixture multinomial models by gradient ascent and greedy search
%PARAMETERS:
%   y - data
%   p1 - first component
%   p2 - second component
%RETURN:
%   alpha - coefficient from gradient ascent
%   besta - coefficient from greedy search
%   sola - coefficient from linear equation solution
%DESCRIPTION:
%
%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   05-MAR-2004 Modified TINGZ
%   04-NOV-2004 Modified TINGZ
%       - add comments
%       - change function name tz_fitmnom --> tz_fitmnom2_ga

zeroc=find((p1+p2)==0);
y(zeroc)=[];
p1(zeroc)=[];
p2(zeroc)=[];

alpha=0.5;
besta=0.5;
mine=1e-5;
maxiter=10000;

adjust=1;
maxadjust=5;

sola=sum(abs(tz_normobjcom(y,0)-p2))/sum(abs(p1-p2));

for i=1:maxiter
    d1=sum(y.*(p1-p2)./(alpha*p1+(1-alpha)*p2));
    d2=-sum(y.*(p1-p2).^2./(alpha*p1+(1-alpha)*p2).^2);
    
    dalpha=d1/d2;
    
    if alpha-dalpha>1
        alpha=(1+alpha)/2;
    else if alpha-dalpha<0
            alpha=alpha/2;
            
        else
            alpha=alpha-dalpha;
        end
    end
    
    if(abs(dalpha)<mine)
        if (alpha>=0 & alpha<=1) | adjust>=maxadjust
            break;
        end
        
        if alpha<0
            alpha=0.01*adjust/maxadjust;
            adjust=adjust+1;
        else 
            alpha=0.99*adjust/maxadjust;
            adjust=adjust+1;
        end
        
    end
    
    
    % %    loglk=sum(y.*log(alpha*p1+(1-alpha)*p2));
    %     
    if i==10000
        warning('max iteration reached');
    end
    
end

prec=100;

zeroy=find(y==0);
tp1=p1;
tp2=p2;
ty=y;
tp1(zeroy)=[];
tp2(zeroy)=[];
ty(zeroy)=[];

for i=0:prec
    a=i/prec;
    if(any(a*tp1+(1-a)*tp2==0))
        loglk(i+1)=-Inf;
    else
        loglk(i+1)=sum(ty.*log((a*tp1+(1-a)*tp2)));
    end
end

[maxlk,besta]=max(loglk);
besta=(besta-1)/100;
%besta=0.5;
if(alpha>1)
    if besta ~=1
        warning('out of range 1');
        besta
    end
    alpha=1;
    break;
end

if(alpha<0)
    if besta ~=0
        warning('out of range 0');
        besta
    end
    alpha=0;
    break
end

