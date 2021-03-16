function regmodel = tz_trainclassif_pslid(x,y,t,regmeth)
%TZ_TRAINCLASSIF_PSLID Obsolete.

%function regmodel = tz_classif_pslid(x,y,t,regmeth)
%OVERVIEW
%   
%PARAMETERS
%   x - 
%   y - 
%   t - 
%   regmeth - 
%       1 - lda
%       2 - svm
%       3 - bpnn      
%RETURN
%   regmodel - 
%DESCRIPTION
%   
%HISTORY
%   22-Jun-2005 Initial write TINGZ
%SEE ALSO
%   

switch regmeth
case 1
    regmodel = tz_ldareg(x,y);
case 2
    if isempty(t)
        clear t;        
        t.norm=1;   
        t.stop=0;   
        t.C_values = 20;    
        t.rbf_levels = 7;  
        t.model_types = maxwin;  
        t.tutor  = smosvctutor;
    end
    y = tz_label2post(y); y(find(y==0)) = -1;
    regmodel = tz_svmreg(x,y,t);    
case 3
    if isempty(t)
        clear t;
        t.norm=2;
        t.epochs=300;
        t.bp=1000;
        t.hidden=20;
        t.stop=1;
    end
    y=tz_label2post(y)*0.8+0.1
    regmodel = tz_bpnnreg(x,y,t);
end

regmodel.postp.ctg=1;
