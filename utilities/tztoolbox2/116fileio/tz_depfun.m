function list = tz_depfun(filename,lim)
%TZ_DEPFUN Check dependencies for a function.
%   LIST = TZ_DEPFUN(FILENAME) returns a list of functions that are
%   required for function FILENAME.
%   
%   LIST = TZ_DEPFUN(FILENAME,LIM) can filter the dependent functions.
%   Only functions with path containing substring LIM are listed.

%   28-Jun-2005 Initial write TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University  

fulllist=depfun(filename);
list=fulllist;
if ~exist('lim','var')
%     list=fulllist;
    return;
end

k=1;
% list={};
for i=1:length(fulllist)
    if isempty(findstr(fulllist{i},lim))
        list(k)=[];
        k=k-1;
    end
    k=k+1;
end