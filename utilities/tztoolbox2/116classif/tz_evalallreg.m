function y = tz_evalallreg(x,regmodel,postp)
%TZ_EVALALLREG Obsolete.

%function y = tz_evalallreg(x,regmodel)
%OVERVIEW
%   
%PARAMETERS
%   x - 
%   regmodel - 
%   postp - 
%RETURN
%   y - 
%DESCRIPTION
%   
%HISTORY
%   12-Jul-2005 Initial write TINGZ
%SEE ALSO
%   

error(tz_genmsg('of','tz_evalallreg','ml_evalreg'));

if isfield(regmodel,'succ')
     if regmodel.succ==0
         warning('An invalid model is taken.');
         y=0;
         return;
     end
end

if nargin>2
    regmodel.postp=postp;
end

switch regmodel.modelname
case 'bpnn'
    y=tz_evalbpnnreg(x,regmodel);
case 'svm'
    y=tz_evalsvmreg(x,regmodel);
case 'lda'
    y=tz_evalldareg(x,regmodel);
end