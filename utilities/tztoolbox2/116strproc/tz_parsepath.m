function dirs=tz_parsepath(fullpath)
%TZ_PARSEPATH Decomposes full path into parts.
%   TZ_PARSEPATH(FULLPATH) returns a cell array of separated directory or
%   file names in FULLPATH. Currently it only works for unix.


%   23-OCT-2004 Inintial write TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University 

if length(fullpath)<=1
    dirs=fullpath;
    return
end

dirs={};

if fullpath(1)=='/'
    dirs{1}='/';
else
    fullpath=['/' fullpath];
end

if fullpath(end)~='/'
    fullpath(end+1)='/';
end

dirpos=find(fullpath=='/');

for i=1:length(dirpos)-1
    dirs{end+1}=fullpath(dirpos(i)+1:dirpos(i+1)-1);
end
