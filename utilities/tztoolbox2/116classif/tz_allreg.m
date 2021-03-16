function regmodel = tz_allreg(x,y,t,regmeth)
%TZ_ALLREG Obsolete

%function regmodel = tz_allreg(x,y,t,regmeth)
%OVERVIEW
%   
%PARAMETERS
%   x - 
%   y - 
%   t - 
%   regmeth - regression method
%       1 - lda
%       2 - svm
%       3 - bpnn
%RETURN
%   regmodel - trained model
%DESCRIPTION
%   
%HISTORY
%   22-Jun-2005 Initial write TINGZ
%SEE ALSO
%   

error(tz_genmsg('of','tz_allreg','ml_regress'));

switch regmeth
case 1
    regmodel = tz_ldareg(x,y);
case 2
    if isempty(t)      
        t.norm=1;   
        t.stop=0;   
        t.C_values = 20;    
        t.rbf_levels = 7;  
        t.model_types = maxwin;  
        t.tutor  = smosvctutor;
    end
    regmodel=tz_svmreg(x,y,t);
case 3
    if isempty(t)
        t.norm=2;
        t.epochs=300;
        t.bp=1000;
        t.hidden=20;
        t.stop=1;
    end
    regmodel = tz_bpnnreg(x,y,t);
end