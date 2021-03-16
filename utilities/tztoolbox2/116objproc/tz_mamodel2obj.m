function obj=tz_mamodel2obj(axln,ds,imgsize)
%TZ_MAMODEL2OBJ Obsolete. Convert object medial axis and width to an object.
%   OBJ = TZ_MAMODEL2OBJ(AXLN,DS,IMGSIZE)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function obj=tz_mamodel2obj(axln,ds,imgsize)
%
%OVERVIEW:
%   convert object medial axis and width to an object
%PARAMETERS:
%   axln - medial axis
%   ds - width
%   imgsize - the size of the image where the object is located
%             empty for automatic image size
%RETURN:
%   obj - object
%DESCRIPTION:
%   
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   05-NOV-2004 Modified TINGZ
%       - add comments
%       - change function name tz_model2img --> tz_mamodel2obj

nholds=length(axln);
axln=axln-min(axln)+1;

if isempty(imgsize)
    axlni=axln;
    dsi=ds;
    imgsize=[max(axln),max(ds)];
else
    xi=(0:imgsize(1)-1)/(imgsize(1)-1)*(nholds-1)+1;
    
    axlni=round(interp1(axln,xi));
    dsi=round(interp1(ds,xi));
    
    axlni=axlni-min(axlni)+1;
end

k=1;
for i=1:length(axlni)
    
    offset=round((imgsize(2)-dsi(i))/2);
    
    for j=1:dsi(i)
        obj(k,:)=[i,axlni(i)+j-1+offset];
        k=k+1;
    end
end

obj(:,3)=1;

%objimg=tz_objimg(obj,[]);
