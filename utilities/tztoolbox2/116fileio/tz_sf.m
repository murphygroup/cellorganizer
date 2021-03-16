function output=tz_sf(str)
%TZ_SF Search files under /home/tingz/matlab
%   LIST = TZ_SF(STR) returns a list of file names matching STR in 
%   /home/tingz/matlab. LIST is a cell array.
%   
%   See also TZ_SEARCHFILE, TZ_SF2, TZ_SFM

%   ??-???-2004 Initial write T. Zhao
%   02-NOV-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

output2=tz_searchfile(str,'/home/tingz/matlab');
if nargout==1
    output = output2;
end

if isempty(output2)
    disp('Nothing is found');
else
    for i=1:length(output2)
        disp(output2{i});
    end
end
