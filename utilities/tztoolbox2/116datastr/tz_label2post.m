function post=tz_label2post(label,N)
%TZ_LABEL2POST Obsolete.
%
%See also ML_LABEL2POST

%function post=tz_label2post(label)
%
%OVERVIEW:
%   convert label to one-of-N coding
%PARAMETERS:
%   label - a vector
%RETURN:
%   post - a matrix
%DESCRIPTION:
%   for the return format of kmeans
%
%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   04-NOV-2004 Modified TINGZ
%       - add comments
%   10-FEB-2005 Modified TINGZ
%       - change algorithm to increase the speed

error(tz_genmsg('of','tz_label2post','ml_label2post'));

if any(label==0)
    error('label could not be zero')
end

if ~isvector(label)
    error('label should be a vector')
end

% tic
% post=zeros(length(label),max(label));
% 
% for i=1:length(label)
%     post(i,label(i))=1;
% end
% toc
maxlabel=max(label);
id=eye(maxlabel);
if nargin==2  
    id=[id,zeros(size(id,1),N-size(id,2))];
end   

post=id(label,:);