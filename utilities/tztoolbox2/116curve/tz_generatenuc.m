function [ds,axln]=tz_generatenuc(lnparas,dsparas,len)
%TZ_GENERATENUC Generate nucleus shape.
%   DS = TZ_GENERATENUC(LNPARAS,DSPARAS,LEN) returns the widths of the
%   gerenrated nucleus along its median axis. LNPARAS contains parameters
%   for generating the median axis and DSPARAS contains parameters for 
%   generating the widths. Both are matrices and each row are parameters
%   for one training nucleus. LEN is a vector of the height of the nucleus.
%   
%   [DS,AXLN] = TZ_GENERATENUC(...) also returns the median axis.
%
%   See also TZ_FITAXLN, TZ_FITDS

%   ??-???-???? Initial write T. Zhao
%   04-NOV-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% sdspara=mvnrnd(mean(dsparas),cov(dsparas),1);
% 
% while sdspara(1)<1
%     sdspara=mvnrnd(mean(dsparas),cov(dsparas),1);
% end
% 
% slnpara=mvnrnd(mean(lnparas),cov(lnparas),1);
% while(any(slnpara([2 4 5 7])<0))
%      slnpara=mvnrnd(mean(lnparas),cov(lnparas),1);
% end

paras=[dsparas,lnparas];
spara=mvnrnd(mean(paras),cov(paras),1);
while spara(1)<1 | any(spara([7 9 10 12])<0)
    spara=mvnrnd(mean(paras),cov(paras),1);
end
sdspara=spara(1:5);
slnpara=spara(6:end);

x=(0:len-1)/(len-1)*2-1;

ds=tz_parabola(sdspara,x);
axln=tz_gaussian(slnpara,x);

plot(1:len,axln-ds/2,'.');
hold on
plot(1:len,axln+ds/2,'.');
 
plot(1:len,axln,'x');
hold off
 
axis('equal')
 
x=1:len;
figure
axln2=axln+ds/2;
plot([x,x(length(x):-1:1),x(1)],[axln-ds/2,axln2(length(x):-1:1),axln(1)-ds(1)/2]);
axis([min(x)-10,max(x)+10,min(axln-ds/2)-10,min(axln-ds/2)+max(x)+10])
axis('equal')