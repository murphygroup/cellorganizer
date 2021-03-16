function fullname = tz_getmodelfullname(modelname)
%TZ_GETMODELFULLNAME Obsolete

%function fullname = tz_getmodelfullname(modelname)
%OVERVIEW
%   
%PARAMETERS
%   modelname - 
%RETURN
%   fullname - 
%DESCRIPTION
%   
%HISTORY
%   25-Jun-2005 Initial write TINGZ
%SEE ALSO
%   

error(tz_genmsg('of','tz_getmodelfullname','ml_getmodelfullname'));

switch(modelname)
case 'bpnn'
    fullname='neural network';
case 'svm'
    fullname='support vector machine';
case 'lda'
    fullname='linear discriminant analysis';
otherwise
    fullname=modelname;
end