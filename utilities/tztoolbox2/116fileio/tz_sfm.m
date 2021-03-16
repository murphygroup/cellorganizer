function tz_sfm(str)
%TZ_SF Search files under ~/matlab
%   LIST = TZ_SFM(STR) returns a list of file names matching STR in 
%   ~/matlab. LIST is a cell array.
%   
%   See also TZ_SEARCHFILE, TZ_SF, TZ_SF2

%   ??-???-2004 Initial write TINGZ
%   02-NOV-2004 Modified TINGZ
%       - add comments

tz_searchfile(str,'~/matlab');