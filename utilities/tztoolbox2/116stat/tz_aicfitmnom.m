function [x,selmodel,aic]=tz_aicfitmnom(y,p)
%TZ_AICFITMNOM Infer multinomial mixture by AIC criteria.
%   X = TZ_AICFITMNOM(Y,P)
%   
%   [X,SELMODEL,AIC] = TZ_AICFITMNOM(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [x,selmodel,aic]=tz_aicfitmnom(y,p)
%OVERVIEW:
%   find best multinomial mixture by aic
%PARAMETERS:
%   y - data
%   p - multinomial prameters mxk
%RETRUN:
%   x - compostion in each component qxk 
%   selmodel - select components 1xq
%   aic - aic value
%DESCRIPTION:
%   this function use tz_redistrmnom to fit multinomial model. It's local maximum.
%
%HISTORY:
%   ??-???-???? Initial write TINGZ

sel=[];
pold=[];
pnew=[];
aic=-Inf;
nmaxmix=size(p,1);
unsel=1:nmaxmix;

k=1;

selmodel=[];
x=y;
while(~isempty(unsel))
    sel=0;
   
    pold=pnew;
    for i=unsel
        pcurr=[pold;p(i,:)];
        if size(pcurr,1)>1
            [tempx,logll]=tz_redistrmnom(y,pcurr);
        else
            tempx=x;
            logll=tz_mnomlogl(x,pcurr);
        end
        nparam=size(pcurr,1)*(size(pcurr,2)-1);
        if aic<logll-nparam
            pnew=pcurr;
            aic=logll-nparam;
            sel=i;
            x=tempx;
        end
        
    end
    
    if sel==0
        break
    else
        selmodel=[selmodel sel];
    end
    unsel(unsel==sel)=[];
end
        
if isempty(selmodel)
    selmodel=1:nmaxmix;
    pcurr=p;
    pnew=p;
    [tempx,logll]=tz_redistrmnom(y,pcurr);
    aic=logll-nmaxmix;
    
    while(length(selmodel)>1)
        unsel=[];
        for i=selmodel
           pcurr=pnew;
           pcurr(selmodel==i,:)=[];
           [tempx,logll]=tz_redistrmnom(y,pcurr);
           nparam=size(pcurr,1)*(size(pcurr,2)-1);
           if aic<logll-nparam
               unsel=i;
               x=tempx;
               aic=logll-nparam;
               
           end
        end
        if isempty(unsel)
            break;
        else
            selmodel(selmodel==unsel)=[];
            pnew(selmodel==unsel,:)=[];
        end
    end
    if(length(selmodel)==1)
        x=y;
    end
end