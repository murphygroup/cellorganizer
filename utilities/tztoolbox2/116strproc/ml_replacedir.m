function fullpath2 = ml_replacedir(fullpath,newdir)
%ML_REPLACEDIR Replace the directory in a full path.
%   FULLPATH2 = ML_REPLACEDIR(FULLPATH,NEWDIR) returns a string that
%   represents the full path of the changed full path FULLPATH. So the
%   new full path becomes NEWDIR/file.

%   10-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

[olddir,file,ext] = fileparts(fullpath);

if(newdir(end) ~= filesep)
    newdir = [newdir filesep];
end

fullpath2 = [newdir file ext];