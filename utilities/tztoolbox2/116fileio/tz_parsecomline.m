function [s,comments] = tz_parsecomline(s,comments)

%TZ_PARSECOMLINE Obsolete

%function [s,comments] = tz_parsecomline(s,comments)
%OVERVIEW
%   parse comments
%PARAMETERS
%   s - 
%   comments - 
%RETURN
%   s - 
%   comments - 
%DESCRIPTION
%   
%HISTORY
%   16-Apr-2005 Initial write TINGZ
%SEE ALSO
%   

error(tz_genmsg('of','tz_parsecomline','tz_parsecmtline'));

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
        s.comparts(i)='1';
        break;
    end
end

comments(1)=[];