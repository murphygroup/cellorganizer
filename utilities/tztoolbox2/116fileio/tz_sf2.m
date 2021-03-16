function list=tz_sf2(str)
%TZ_SF2 Search files under current directory.
%   LIST = TZ_SF2(STR) returns a list of file names matching STR in 
%   current directory. LIST is a cell array.
%   
%   See also TZ_SEARCHFILE, TZ_SF, TZ_SFM

%   ??-???-2004 Initial write TINGZ
%   02-NOV-2004 Modified TINGZ
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

list=tz_searchfile(str,'./');