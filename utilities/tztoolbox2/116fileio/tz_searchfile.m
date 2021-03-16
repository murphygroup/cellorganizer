function list=tz_searchfile(str,sdir,varargin)
%TZ_SEARCHFILE Search files under a directory.
%   TZ_SEARCHFILE(STR,SDIR) returns a list of found file names matching
%   pattern STR under SDIR and its subdirectories. It supports wild cards.
%   
%   TZ_SEARCHFILE(STR,SDIR,OPTIONAL) specifies find options in OPTINAL.
%   
%   See also TZ_SF, TZ_SF2, TZ_SFM

%   ??-???-2004 Initial write TINGZ
%   02-NOV-2004 Modified TINGZ
%       - add comments
%       - add parameter sdir
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('2 or 3 arguments are required')
end

pstar=find(str=='*');
str2=str;
if ~isempty(pstar)
    str2=str(1:pstar(1)-1);
    for i=1:length(pstar)-1
        str2=[str2 '\' str(pstar(i):pstar(i+1)-1)];
    end
    str2=[str2,'\',str(pstar(end):end)];
end

cmd = ['find ' sdir ' ' ' -name ' str2 ' -print'];
if nargin==3
    cmd = [cmd ' ' varargin{1}];
end
[status,output] = unix(cmd);

if isempty(output)
    list=[];
    return;
end

linepos=[0 find(output==10)];
   
k=1;
for i=1:length(linepos)-1
    filepath=output(linepos(i)+1:linepos(i+1));
    if isempty(findstr(filepath,'Permission denied')) & ~isempty(filepath)
        list{k,1}=filepath(1:end-1);
        k=k+1;
    end
end