function [p_value,both] = tz_multest2(f1,f2,Method,TestMethod,dependence,sigvalue,pvalues)
%TZ_MULTEST2 Obsolete.
%
%See also ML_MULTEST2

%[p_value,both] = tz_multest2(f1,f2,Method,TestMethod,dependence,sigvalue,pvalues)
%OVERVIEW:
%   Multiple testing
%PARAMETERS:
%   Method - mutltiple testing method
%   Testmethod - test method
%   f1 - sample 1
%   f2 - sample 2
%   dependence - dependent or not, for BH method
%   sigvalue - significant level
%   pvalues - calculated pvalue
%RETURN:
%   both - A vector indicating which features are changed
%   p_value - pvalues
%DESCRITPION:
%   
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   14-JUL-2003 Modified TINGZ
%   09-NOV-2003 Modified TINGZ 
%   14-DEC-2003 Modified TINGZ
%   26-JUL-2004 Modified TINGZ
%       - change the order of output
%       - change the order of parameters
%   03-AUG-2004 Modified TINGZ
%       - Add p-values to the last row of the output both

error(tz_genmsg('of','tz_multest2','ml_multest2'));

feature_size=size(f1,2);

if dependence==1
    Cm=sum(1./[1:feature_size]);
else
    Cm=1;
end

switch(TestMethod)
case 'ready'
    p=pvalues;
    
case 'ttest2'
    for k=1:feature_size
        s1=f1(:,k);
        s2=f2(:,k);
        if isempty(tz_constidx([s1;s2]))
            [h, p(k), ci]=eval(strcat(TestMethod,'(s1,s2)'));
        else
            p(k)=1;
        end
    end
otherwise
    for k=1:feature_size
        s1=f1(:,k);
        s2=f2(:,k);
        
        if isempty(tz_constidx([s1;s2]))
            p(k)=eval(strcat(TestMethod,'(s1,s2)'));
        else
            p(k)=1;
        end
    end
end

[sorted_p index]=sort(p);
rejections=zeros(1,feature_size);

switch(Method)
case 'Un'
    p_value=sorted_p(1);
    both=[index;sorted_p<=sigvalue;sorted_p];
case 'BH'

    for(i=feature_size:(-1):1)
        if(sorted_p(i)<=(i*sigvalue/feature_size/Cm))
            rejections(1:i)=1;
            break;
        end
    end
    
    p_value=min(1,min(sorted_p*Cm*feature_size./(1:feature_size)));
    both=[index; rejections];
case 'Bn'
    t=sigvalue/feature_size;
    p_value=min(sorted_p(1)*feature_size,1);
    both=[index; sorted_p<=t; sorted_p];
case 'Hm'
    for(i=1:feature_size)
        if(sorted_p(i)>sigvalue/(feature_size-i+1))
            rejections=rejections+1;
            rejections(i:feature_size)=0;
            break;
        end
    end
    
    p_value=min(sorted_p(1)*feature_size,1)
    both=[index; rejections; sorted_p]
end