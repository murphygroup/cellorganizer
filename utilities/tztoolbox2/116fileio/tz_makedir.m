function tz_makedir( dirname,os)

%TZ_MAKEDIR Obsolete

%function tz_makedir(dirname,os)
%
%OVERVIEW:
%   make a directory under current directory
%PARAMETERS:
%   dirname - the new directory name
%   os - operating system
%RETURN:
%
%DESCRIPTION:
%   You'd better use mkdir function
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   21-JUL-2003 Modified TINGZ
%   12-JAN-2004 Modified TINGZ

warning(tz_genmsg('of','tz_makedir'));

if(~exist('os','var'))
    os='unix'
end

status=0;

switch(os)
case 'unix',
    if( ~exist( dirname, 'dir'))
        command = ['mkdir ' dirname]
        status = unix( command);
    end
case 'dos',
    if( ~exist( dirname, 'dir'))
        command = ['mkdir ' tz_dirunix2dos(dirname)]
        status = dos(command);
    end
end

if( status ~= 0)
    error(['Cannot create directory: ' dirname]);
end
