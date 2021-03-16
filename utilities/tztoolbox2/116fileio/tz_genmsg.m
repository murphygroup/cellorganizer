function msg = tz_genmsg(msgtype,varargin)
%TZ_GENMSG Generates message.
%   MSG = TZ_GENMSG(MSGTYPE) automatically generates message with
%   type MSGTYPE:
%       'of' - obsolete function
%       'ip' - parameter initialization
%   MSG = TZ_GENMSG('cv',VARNAME,FUNNAME) generates a message for a 
%   nonexisted variable with the name VARNAME, which is supposed to be
%   generated from the script FUNNAME.

%   18-Apr-2005 Initial write TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University  

switch msgtype
    case 'of' %obsolete function
        if nargin >= 3
            replaceMsg = ['Please use ' varargin{2} ' instead.'];
        else
            replaceMsg = 'Currently there is no replacement.';
        end

        msg=sprintf('%s\n%s',...
            ['The function ' varargin{1} ' is out of date.'],...
            replaceMsg);
    case 'ip' %parameter initialization
        msg=sprintf('%s\n%s\n%s',...
            ['if ~exist(''param'',''var'')'], ...
            ['    param = struct([])'], ...
            ['end']);
    case 'cv' %check variable
        if nargin < 3
            error('Exactly 3 arguments are required')
        end
        
        varname = varargin{1};
        funname = varargin{2};
        msg = sprintf('%s\n%s',...
            ['The variable ' varname ' does not exist.'], ...
            ['Please run ' funname ' first or check ' funname '.']);
    case 'bg' %show testing message
        msg = sprintf('%s',...
            ['Possible bug in ' varargin{1} ]);
    otherwise
        error(['Unrecognized message type:' msgtype]);
end
