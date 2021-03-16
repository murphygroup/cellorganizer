function features=tz_combfeats2mcf(combfeats,combclass)
%TZ_COMBFEATS2MCF Obsolete
%
%See also ML_COMBFEATS2MCF

%function features=tz_combfeats2mcf(combfeats,combclass)
%   
%OVERVIEW:
%   convert combined features to mcf
%PRAMETERS:
%   combfeats - combined feature matrix
%   combclass - combined class label
%RETURN
%   features - mcf
%DESCRIPTION
%
%HISTORY
%   ??-???-???? Initial write TINGZ
%   31-OCT-2004 Modified TINGZ

error(tz_genmsg('of','tz_combfeats2mcf','ml_combfeats2mcf'));

features=tz_findclass([combfeats,combclass]);
nclass=length(features);

for i=1:nclass
    features{i}=features{i}(:,1:end-2);
end