function tz_setenv(varargin)
%TZ_SETENV Set enviromental variables.
%   TZ_SETENV('dbg',1) turns debug status on and TZ_SETENV('dbg',0) turns
%   it off.
%   TZ-SETENV('upd',1) turns updating status on and TZ_SETEVN('upd',1)
%   turns it off.
%   TZ_SETENV('rid',VALUE) sets result id to value VALUE.
%   This funciton also supports set multiple variables at the same time.
%   For example:
%       TZ_SETENV('dbg',1,'upd',0)
%
%   See also TZ_GETENV

%   05-DEC-2004 Initial write  Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU


if mod(nargin,2)~=0
    error('input missed');
end

nvar=nargin/2;
markfile='tz_config.m';
envdir=which(markfile);
envdir=envdir(1:end-length(markfile)-1);
envfile=[envdir '/' 'env.mat'];

if exist(envfile,'file')
    load(envfile);
end

for i=1:nvar
    option=varargin{2*i-1};
    switch option
    case 'dbg'
        tzv_debug=varargin{2*i};
    case 'upd'
        tzv_update=varargin{2*i};
    case 'rid'
        tzv_resultid=varargin{2*i};
    otherwise
        warning(['wrong option: ' option]);
    end
end

save(envfile,'tzv_*');
        
        


