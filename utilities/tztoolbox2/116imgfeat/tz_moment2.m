function mom = tz_moment2(img)
%TZ_MOMENT2 Obsolete.

%function mom = tz_moment2(img)
%OVERVIEW   
%   Calculate central moments of an image. Same as tz_moment
%PARAMETERS
%   img - input image
%RETURN
%   mom - moments
%DESCRIPTION
%   
%HISTORY
%   24-Mar-2005 Initial write TINGZ
%SEE ALSO
%   tz_bwmoment, tz_moment

error('Function tz_moment2 is out of date. Please use ml_moment2');

if sum(img(:))==0
    mom=[];
    warning('black image');
    return;
end

[x,y]=find(img>0);
weights=img(find(img>0));
S=sum(weights);
centroid=sum([x.*weights,y.*weights],1)/S;
cpos=tz_addrow([x,y],-centroid);
% cpos=[x,y];

mu00=1;
mu11=sum(weights.*cpos(:,1).*cpos(:,2))/S;
mu20=sum(weights.*cpos(:,1).^2)/S;
mu02=sum(weights.*cpos(:,2).^2)/S;
% mu30=sum(weigths.*cpos(:,1).^3)/S;
% mu03=sum(weigths.*cpos(:,2).^3)/S;
mom = struct('cx',centroid(1),'cy',centroid(2),'mu00',mu00,'mu11',mu11,'mu20', mu20,'mu02',mu02);
