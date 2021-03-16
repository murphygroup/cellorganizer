function [features,names]=tz_loadobj_sd(objdir)

%TZ_LOADOBJ_SD Obsolete

%function [features,names]=tz_loadobj_sd(objdir)
%OVERVIEW:
%   Load obj level data from a directory
%PARAMETERS:
%   objdir - directory
%RETURN:
%   features - object level features
%   names - group names
%DESCRIPTION:
%   no idea about the use of m
%
%HISTORY:
%   11-FEB-2004 Initial write TINGZ
%   12-FEB-2004 Modified TINGZ

dirlist=xc_readdir(objdir);

file_no=length(dirlist);
m=8;

names{1}=dirlist{1}(1:m);
features{1}=[];
k=1

for i=1:file_no
    if all(names{k}==dirlist{i}(1:m))
        load([objdir '/' dirlist{i}]);
        features{k}=[features{k};feats];
    else
        k=k+1
        names{k}=dirlist{i}(1:m);
        features{k}=[];
    end
end