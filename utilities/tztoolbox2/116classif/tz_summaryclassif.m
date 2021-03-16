function [ncm,pcm,avgacc]=tz_summaryclassif(tlabel,elabel)
%TZ_SUMMARYCLASSIF Obsolete

%function [ncm,pcm,avgacc]=tz_summaryclassif(tlabel,elabel)
%
%OVERVIEW:
%   Summarize classification results
%PARAMETERS:
%   tlabel - true label mx1
%   elabel - estimated label nx1
%RETRUN:
%   ncm - confusion matrix by numbers
%   pcm - confusion matrix by percentage
%   avgacc - average accuracy
%DESCRIPTION:
%
%HSITORY:
%   17-MAY-2004 Modified TINGZ
%       - return confusion matrix of numbers
%SEE ALSO
%

error(tz_genmsg('of','tz_summaryclassif','ml_summaryclassif'));

%find distinguished classes
caclass=tz_findclass(tlabel);

nclass=length(caclass);
nsample=length(elabel);

ncm=zeros(nclass,nclass);

for i=1:nsample
    ncm(tlabel(i),elabel(i))=ncm(tlabel(i),elabel(i))+1;
end

pcm=ncm./repmat(sum(ncm,2),1,nclass)*100;
avgacc=sum(tlabel==elabel)/nsample;