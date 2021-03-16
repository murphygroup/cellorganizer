function [alpha,selmodel,aic]=tz_aicfitmnomk(y,p)
%TZ_AICFITMNOMK Infer multinomial mixture by AIC criteria.
%   ALPHA = TZ_AICFITMNOMK(Y,P)
%   
%   [ALPHA,SELMODEL,AIC] = TZ_AICFITMNOMK(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [alpha,selmodel,aic]=tz_aicfitmnomk(y,p)
%OVERVIEW:
%   find best multinomial mixture by aic
%PARAMETERS:
%   y - data
%   p - multinomial prameters mxk
%RETRUN:
%   alpha - compostion in each component 1xq 
%   selmodel - select components 1xq
%   aic - aic value
%DESCRIPTION:
%   this function use tz_fitmnomk to fit multinomial model.
%   the likelihood is global maxima
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


while(~isempty(unsel))
    pold=pnew;
    sel=0;
    for i=unsel
        pcurr=[pold;p(i,:)];
        nparam=size(pcurr,1)*(size(pcurr,2)-1);
        [tmpalpha,loglk]=tz_fitmnomk(y,pcurr);
        tmpaic=loglk-nparam;
        if tmpaic>aic
            alpha=tmpalpha;
            pnew=pcurr;
            aic=tmpaic;
            sel=i;
        end
    end
    if sel==0
        break;
    else
        selmodel=[selmodel,sel];
        unsel(unsel==sel)=[];
    end
end