function s = tz_genargcheck(ninput,noptional)
%TZ_GENARGCHECK Genearte codes for argument checking
%   S = TZ_GENARGCHECK(NINPUT) returns a string used for checking number
%   of arguments to see if it is matched with NINPUT.
%   
%   S = TZ_GENARGCHECK(NINPUT,NOPTIONAL) also considers number of 
%   optional arguments
%  
%   See also

%   12-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required')
end

s='';

if ninput == 0
    return;
end

s = tz_addstrline(s,['if nargin < ' num2str(ninput)]);
matlabIndent = '    ';

switch noptional     
case 0
    s = tz_addstrline(s,[matlabIndent,'error(''Exactly ', ...
            num2str(ninput),' argument',pluralstring(ninput), ...
            ' required'');']);
case 1
    s = tz_addstrline(s,[matlabIndent,'error(''',num2str(ninput), ...
            ' or ', num2str(ninput+1),' argument',pluralstring(2), ...
            ' required'');']); 
otherwise
    s = tz_addstrline(s,[matlabIndent,'error(''','At least ', ...
            num2str(ninput),' argument',pluralstring(2), ...
            ' required'');']);
        
end

s = tz_addstrline(s,'end');

function s = pluralstring(n)

if n == 1
    s = ' is';
else
    s = 's are';
end
