function y = tz_testclassif_pslid(sample,regmodel)
%TZ_TESTCLASSIF_PSLID Obsolete.

%function y = tz_testclassif_pslid(sample,regmodel)
%OVERVIEW
%   
%PARAMETERS
%   sample - testing data
%   regmodel - trained classifier
%RETURN
%   y - 
%DESCRIPTION
%   
%HISTORY
%   22-Jun-2005 Initial write TINGZ
%SEE ALSO
%   tz_trainclassif_pslid

switch regmodel.modelname
case 'bpnn'
    y=tz_evalbpnnreg(sample,regmodel);
case 'svm'
    
case 'lda'
    y=tz_evalldareg(sample,regmodel);
end