function s = tz_parsecomments(filename)
%TZ_PARSECOMMENTS Parse comments in a file.
%   S = TZ_PARSECOMMENTS(FILENAME) returns a structure S, which has two
%   string fields:
%       funtitle: 'xxx' (input functionname output)
%       cmtparts: 'xxxxxx' keywords
%   Here x is 0 or 1
%
%   See also TZ_PARSECMTLINE

%   13-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University
  

sfun=tz_parsemfun(filename);

%extract comments
comments=help(filename);
comments = strread(comments,'%s','delimiter','\n','whitespace','');

blanklines=[];
for i=1:length(comments)
    if isempty(tz_strdel(comments{i},' '))
        blanklines=[blanklines,i];
    end
end

if ~isempty(blanklines)
    comments(blanklines)=[];
end

[scomfun,n]=tz_parsemfunlines(comments);

s.funtitle='111';
s.cmtparts='000000';

if ~isempty(sfun.input) | ~isempty(scomfun.input)
    if ~strcmp(tz_cell2str(sfun.input,','),tz_cell2str(scomfun.input,','))
        s.funtitle(1)='0';
    end
end

if ~strcmp(sfun.funname,scomfun.funname)
    s.funtitle(2)='0';
end

if ~isempty(sfun.output) | ~isempty(scomfun.output)
    if ~strcmp(tz_cell2str(sfun.output,','),tz_cell2str(scomfun.output,','))
        s.funtitle(3)='0';
    end
end

if length(comments)==n
    return;
end

comments=comments(n+1:end);

while ~isempty(comments)
    [s,comments]=tz_parsecmtline(s,comments);
end
