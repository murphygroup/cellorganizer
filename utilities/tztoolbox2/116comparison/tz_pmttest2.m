function pvalue = tz_pmttest2(s1,s2)
%TZ_PMTTEST2 Obsolete.
%
%See also ML_PMTTEST2

%function pvalue = tz_pmttest2(s1,s2)
%OVERVIEW:
%   permutation test
%PARAMETERS:
%   s1 - sample 1
%   s2 - sample 2
%RETURN:
%   pvalue - p-value
%DESCRIPTION:
%   50000 trials are made
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   14-JUL-2003 Modified TINGZ 
%   15-DEC-2003 Modified TINGZ

error(tz_genmsg('of','tz_pmttest2','ml_pmttest2'));

t1=reshape(s1,1,size(s1,1)*size(s1,2));
t2=reshape(s2,1,size(s2,1)*size(s2,2));

tobs=abs(mean(t1)-mean(t2));
merges=[t1 t2];
L1=length(t1);
L2=length(t2);
Lm=length(merges);

B=50000;
sum=0;

for i = 1:B
    randorder=tz_randorder(Lm);
    %if any(randorder(1:L1)>L1)
    ro1=randorder([1:L1]);
    ro2=randorder([L1+1:Lm]);
    per1=merges(ro1);
    per2=merges(ro2);
    T=abs(mean(per1)-mean(per2));
    sum=sum+(T>tobs);
    %end
    if  mod(i,10000)==0
        i
    end
end

pvalue=sum/B;
    
