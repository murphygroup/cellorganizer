function img2 = tz_normalize(img)
%TZ_NORMALIZE Unknown.

%function img2 = tz_normalize(img)
%OVERVIEW
%   
%PARAMETERS
%   img - 
%RETURN
%   img2 - 
%DESCRIPTION
%   
%HISTORY
%   28-Jun-2005 Initial write TINGZ
%SEE ALSO
%   

imgsize=size(img);
tcenter=round(imgsize/2);

mom = tz_moment2(img);

