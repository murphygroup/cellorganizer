function tz_updatepslid(filename)
%TZ_UPDATEPSLID Update PSLID matlabshared function.
%   TZ_UPDATEPSLID(FILENAME) copies function with the name FILENAME to
%   PSLID matlab shared directory.

%   20-Jul-2005 Initial write T. Zhao
%   03-Feb-2006 Modified T. Zhao
%       - change pslid server to pslid2
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argumentis required')
end

fullpath=which(filename);
if isempty(fullpath)
    cmd=[];
    comments='no such file, update fialed';
else
    cmd=['scp ' fullpath ... 
            ' tingz@pslid2.cbi.cmu.edu:' ...
            '/usr/local/jakarta-tomcat-5.0.28/webapps/develop/matlabshared'];
    comments=[fullpath ' uploaded to PSLID'];
end

if ~isempty(cmd)
    status=unix(cmd);
    
    if status~=0
        comments=['update fialed'];
    end
    disp(comments);
end
