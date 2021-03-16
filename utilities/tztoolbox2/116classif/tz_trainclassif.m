function regmodel = tz_trainclassif(x,y,t,regmeth)
%TZ_TRAINCLASSIF Obsolete.

%function regmodel = tz_trainclassif(x,y,t,regmeth)
%OVERVIEW
%   
%PARAMETERS
%   x - nxm
%   y - nx1
%   t - structure
%   regmeth - number
%       1 - lda
%       2 - svm
%       3 - bpnn
%RETURN
%   regmodel - 
%DESCRIPTION
%   
%HISTORY
%   12-Jul-2005 Initial write TINGZ
%SEE ALSO
%   tz_allreg

error(tz_genmsg('of','tz_trainclassif','ml_trainclassif'));

switch regmeth
case 1
    regmodel = tz_ldareg(x,y);
case 2
    y = tz_label2post(y); 
    y(find(y==0)) = -1;  
case 3
    y=tz_label2post(y)*0.8+0.1 
end

regmodel = tz_allreg(x,y,t,regmeth);
regmodel.postp.ctg=1;
