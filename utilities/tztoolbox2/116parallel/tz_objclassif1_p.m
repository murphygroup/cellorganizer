function tz_objclassif1_p(k,pfold)

%function tz_objclassif1_p(k)
%OVERVIEW
%   parallel process for object level classification
%PARAMETERS
%   k - number of clusters
%RETURN
%
%DESCRIPTION
%   
%HISTORY
%   05-Apr-2005 Initial write TINGZ
%SEE ALSO
%   

if nargin<2
    pfold=[];
end

tz_pinitpath
load pp_combobjfeats.mat


tic
workdir=[mtpath '/' datadir '/paper2/perm1'] 

[avgcm,avgacc,aic]=tz_objclassif(combobj(:,[1 3 5:end]),combcellidx,combclass,[],{},combobj(:,4),[],...
    10,pfold,'kmeans',k,[],'dist',{},'bpnn',0,{1},'#!objnum',workdir,6,...
    {'35randseeds'})
toc