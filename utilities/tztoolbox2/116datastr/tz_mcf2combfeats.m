function [combfeats,combclass,combcellidx]=tz_mcf2combfeats(features)
%TZ_MCF2COMBFEATS Obsolete.
%
%See also ML_MCF2COMBFEATS

%function [combfeats,combclass,combcellidx]=tz_combfeats_mcf(features)
%
%OVERVIEW:
%   convert mcf to combined features
%PARAMETERS:
%   features - mcf
%RETURN
%   combfeats - combined feature matrix
%   combclass - combined class label
%   combcellidx - combined cell index
%DESCRIPTION:
%
%HISTORY
%   AUG-07-2004 Initial write TING
%   31-OCT-2004 Modified TINGZ
%       - add comments
%       - change function name tz_combinefeats_mcf --> tz_mcf2combfeats

error(tz_genmsg('of','tz_mcf2combfeats','ml_mcf2combfeats'));

nclass=length(features);
combfeats=[];
combclass=[];
combcellidx=[];

for i=1:nclass
    combfeats=[combfeats;features{i}];
    ncell=size(features{i},1);
    combclass=[combclass;zeros(ncell,1)+i];
    combcellidx=[combcellidx;(1:ncell)'];
end
