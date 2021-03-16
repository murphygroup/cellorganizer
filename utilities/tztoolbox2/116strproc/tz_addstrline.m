function str = tz_addstrline(str,strline)
%TZ_ADDSTRLINE adds a line to a string
%   STR = TZ_ADDSTRLINE(STR) adds a line break at the end of STR
%   
%   STR = TZ_ADDSTRLINE(STR,STRLINE) adds STRLINE as a new line of STR
%   
%   See also

%   05-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin<2
    strline = [];
end

if isempty(strline)
    str = sprintf('%s\n',str); %add a line break
else
    if isempty(str)
        str = sprintf('%s',strline);
    else
        str = sprintf('%s\n%s',str,strline);
    end
end