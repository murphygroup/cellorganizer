function featsets=tz_parsefeatset(featset)
%TZ_PARSEFEATSET  Parse feature set string
%   TZ_PARSEFEATSET(FEATSET) separates the feature set string FEATSET into 
%   tokens and returns them in a cell array. Only token following '#'
%   will be taken as a feature set name.
%   
%   Example:
%       tz_parsefeatset('#feat1#feat2') returns {'feat1','feat2'}
%       tz_parsefeatset('feat1#feat2')  returns {'feat2'}
%
%   See also TZ_STRTOK

%   ??-???-???? Initial write TINGZ
%   05-NOV-2004 Modified TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin ~= 1
    error('Exactly 1 argmument is required');
end

sepidx=find(featset=='#');

if isempty(sepidx)
    error('wrong feature set');
end

featsets={};

for i=1:length(sepidx)-1
    featsets{i}=featset(sepidx(i)+1:sepidx(i+1)-1);
end

featsets{end+1}=featset(sepidx(end)+1:end);