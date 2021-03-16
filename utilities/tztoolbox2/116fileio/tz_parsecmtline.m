function [s,comments] = tz_parsecmtline(s,comments)
%TZ_PARSECMTLINE Parse 'tz' style comments.
%   [S2,COMMENTS2] = TZ_PARSECMTLINE(S,COMMENTS) analyzes cell array of 
%   strings, COMMENTS and take the first keyword out. The corresponding
%   postion of S.comparts is set to '1' and returns S2. The the strings 
%   related to the extracted keyword are removed from comments, which 
%   generates the remained string COMMENTS2.
%   
%   See also TZ_PARSECOMMENTS

%   16-Apr-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

keywords=tz_getcmtkeywords;

if isempty(comments)
    s=s;
    comments=comments;
    return;
end

compos=[];
% comments{1}
for i=1:length(keywords)
    
    comline=tz_strdel(upper(comments{1}),' :');
    
    if strcmp(comline,keywords{i})
        s.cmtparts(i)='1';
        break;
    end
end

comments(1)=[];