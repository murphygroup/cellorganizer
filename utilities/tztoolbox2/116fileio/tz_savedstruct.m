function tz_savedstruct(filename,s)
%TZ_SAVEDSTRUCT Obsolete. See ML_SAVESTRUCT.
%   TZ_SAVEDSTRUCT(FILENAME,S) save all fiels in S into the file with
%   file name FILENAME.

%   07-Apr-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_savedstruct','ml_savestruct'));

if nargin < 1
    error('Exactly 1 argumentis required')
end

c=struct2cell(s);
f=fieldnames(s);

savecmd=['save ' filename];

for i=1:length(c)
    eval([f{i} '=' 's.' f{i} ';']);
    savecmd=[savecmd ' ' f{i}];
end

eval(savecmd);

