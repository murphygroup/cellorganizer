function all_features2=tz_selfeat_mcf(all_features,featset,er)
%TZ_SELFEAT_MCF Under construction.

%function all_features2=tz_selfeat_mcf(all_features,featset)
%
%OVERVIEW:
%   select feat set
%PARAMETERS:
%   all_features - mcf
%   featset - see er_selfeat
%RETURN:
%   all_features2 - return mcf
%DESCRIPTION:
%   
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   03-NOV-2004 Modified TINGZ
all_features2=all_features;
nsample=length(all_features(:));
if er==0
    for i=1:nsample
        all_features2{i}=double(all_features{i}(:,featset));
    end
else
    for i=1:nsample
        all_features2{i}=double(er_selfeat(all_features{i},featset));
    end
end