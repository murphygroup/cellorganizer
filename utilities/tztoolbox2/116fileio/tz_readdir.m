function dirname = xc_readdir(expression, tmpdir, os)
%TZ_READDIR Obsolete

% FUNCTION DIR = XC_READDIR(EXPRESSION, TMPDIR)
% same as 'ls expression' in linux.  Support remote server
% @param tmpdir - the tmp directory used

%Last modified by tingz on Jul. 21, 2003
%Last modified by tingz on Jan. 11, 2004
%Last modified by tingz on Apr. 2, 2004
%   Detect os automatically

if(~exist('os','var'))
    if strcmp(computer,'PCWIN')
        os='windows';
    else
        os='unix';
    end
end

switch(os)
case 'unix',
    if (~exist('tmpdir', 'var'))
        tmpdir = '/tmp';
    end
    
    unix(['ls ' expression ' >' tmpdir '/list.xc']);
    dirname = textread([tmpdir '/list.xc'], '%s', 'delimiter', '\n');
    delete([tmpdir '/list.xc']);
case 'windows',
     if (~exist('tmpdir', 'var'))
         tmpdir = 'c:/windows/temp';
     end
%     command = ['dir ' '/b ' expression ' >' tmpdir '\list.xc'];
%     dos(command);
    files=dir(expression);
    for i=1:length(files)
        dirname{i}=files(i).name;
    end
end

%path(path, '/home/shared/matlab');
%pos = findstr(expression, ':');

%save([tmpdir '/list.xc' ],'dirlist','-ASCII');
  
dirname = sort(dirname);
%dirname = {.' dirname{1:end}};