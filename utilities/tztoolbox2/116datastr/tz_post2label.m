function labels=tz_post2label(post)
%TZ_POST2LABEL Obsolete.
%
%See also ML_POST2LABEL

%function labels=tz_post2label(post)
%
%OVERVIEW:
%   convert post to label
%PARAMETERS:
%   post - a matrix
%RETURN:
%   label - a vector
%DESCRIPTION:
%   for the return format of kmeans
%
%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   05-NOV-2004 Modified TINGZ
%       - add comments

error(tz_genmsg('of','tz_post2label','ml_post2label'));

if any(sum(post,2)~=1)
    error('wrong post')
end

[m,labels]=max(post');
labels=labels';