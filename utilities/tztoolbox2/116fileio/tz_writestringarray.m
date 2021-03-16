function tz_writestringarray(outputfile,s)
%TZ_WRITESTRINGARRAY Obsolete

%function tz_writestringarray(outputfile,s)
%OVERVIEW
%   
%PARAMETERS
%   outputfile - 
%   s - 
%RETURN
%
%DESCRIPTION
%   
%HISTORY
%   25-Jun-2005 Initial write TINGZ
%SEE ALSO
%   

error(tz_genmsg('of','tz_writestringarray','ml_writestringarray'));

fid=fopen(outputfile,'w');
for i=1:length(s)
    fprintf(fid,'%s\n',s{i});
end
fclose(fid);