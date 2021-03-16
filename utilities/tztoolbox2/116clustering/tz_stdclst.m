function post2=tz_stdclst(post)
%TZ_STDCLST Obsolete. 
%
%See also ML_STDCLST

%function post2=tz_stdclst(post)
%OVERVIEW:
%   Remove empty clusters
%PARAMETERS:
%   post - original post
%RETURN:
%   post2 - new post
%DESCRIPTION:
%   Clusters from the first classifier may have some empty groups
%
%HISTORY:
%   ??-???-???? Initial write TINGZ

error(tz_genmsg('of','tz_stdclst','ml_stdclst'));

clstsize=sum(post,1);

post2=post;

post2(:,clstsize==0)=[];
