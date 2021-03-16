function cnames=tz_cleandirs(names,cleans)
%TZ_CLEANDIRS Remove strings starting with '.'
%   TZ_CLEANDIRS(NAMES) removes strings starting with '.' in the 
%   string array NAMES and return the new string array.

%   ??-???-???? Initial write TINGZ
%   31-OCT-2004 Modified TINGZ
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin<1
    error('1 or 2 arguments are required.');
end

if ~iscell(names)
    error('The input must be a cell array.');
end

rmi=[];
for i=1:length(names)
    if names{i}(1)=='.'
        rmi=[rmi,i];
    else
        if exist('cleans','var')
            if strcmp(names{i},cleans)
                rmi = [rmi,i];
            end
        end
    end
end

cnames=names;
cnames(rmi)=[];
