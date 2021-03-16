function both = tz_compfeats_u(feature_matrix1,feature_matrix2,method,issorted)
%TZ_COMPFEATS_U Obsolete.
%
%See also ML_COMPFEATS_U

%function both = tz_compfeat_u(method,feature_matrix1,feature_matrix2,issorted)
%OVERVIEW:
%   univariate test
%PARAMETERS:
%   feature_matrix1 - the first matrix n1xp
%   feature_matrix2 - the seconde matrix n2xp
%   method - test method
%   issorted - sort the pvalues or not; 0 or 1
%RETURN:
%   both - pvalues 2xp
%
%HISTORY
%   ??-???-???? Initial write TINGZ
%   09-JUN-2003 Modified TINGZ
%   14-JUN-2003 Modified TINGZ
%   14-DEC-2003 Modified TINGZ
%   18-FEB-2004 Modified TINGZ

error(tz_genmsg('of','tz_compfeats_u', 'ml_compfeats_u');

if ~exist('issorted','var')
    issorted=1;
end

both=1:size(feature_matrix1,2);
pvalues=[];

for featindex=1:size(feature_matrix1,2)
    s1=feature_matrix1(:,featindex);
    s2=feature_matrix2(:,featindex);
    if isempty(tz_constidx([s1;s2]))
        pvalue=eval(strcat(method,'(s1,s2)'));
    else
        pvalue=1;
    end
    pvalues=[pvalues pvalue];
        
end

both=[both;pvalues];

if issorted==1
    conf = both(2,:);
    [sorted_conf index] = sort(conf');
    both = [index sorted_conf]';
end

