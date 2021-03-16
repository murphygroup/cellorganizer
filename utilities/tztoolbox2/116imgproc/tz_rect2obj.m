function objimg=tz_rect2obj(rcimg,obj)
%TZ_RECT2OBJ Morph a rectangle to an object.
%   OBJIMG = TZ_RECT2OBJ(RCIMG,OBJ)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_rect2obj(img,obj)
%
%OVERVIEW:
%   morph an image to an object
%PARAMETERS:
%   rcimg - input image
%   obj - object
%RETURN:
%   objimg - output image
%DESCRIPTION:
%   
%HISTORY:
%   ??-???-???? Intial write TINGZ
%   05-NOV-2004 Modified TINGZ
%       -add commetns

objrange=tz_findrange(obj);

obj(:,2)=obj(:,2)-min(obj(:,2))+1;

rangey=[min(obj(:,1)),max(obj(:,1))];

rcimg=imresize(rcimg,[rangey(2)-rangey(1)+1,size(rcimg,2)],'bilinear');

rcimg=[zeros(1,size(rcimg,2));rcimg;zeros(1,size(rcimg,2))];

objimg=zeros(objrange);
k=1;
for i=rangey(1):rangey(2)
    hln=obj(obj(:,1)==i,2);
    if ~isempty(hln)
        subimg=rcimg(k:k+2,:);
        subimg=imresize(subimg,[size(subimg,1),length(hln)]);
        objimg(k,hln)=subimg(2,:);
%         ns=(1:length(hln))/length(hln)*size(rcimg,2);
%         objimg(k,hln)=interp1(rcimg(k,:),ns);
    end
    k=k+1;
end
