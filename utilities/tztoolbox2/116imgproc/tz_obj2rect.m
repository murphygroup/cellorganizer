function [img,objimg]=tz_obj2rect(obj,imgsize)
%TZ_OBJ2RECT Morph an objct to an rectangle.
%   IMG = TZ_OBJ2RECT(OBJ,WNDSIZE)
%   
%   [IMG,OBJIMG] = TZ_OBJ2RECT(...)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function img=tz_obj2rect(obj,size)
%
%OVERVIEW:
%   morph an objet to an image
%PARAMETERS:
%   obj - input object
%   imgsize - size of the image
%RETURN:
%   img - output image
%   objimg - original object image
%DESCRIPTION:
%   
%HISTORY:
%   ??-???-???? Intial write TINGZ
%   05-NOV-2004 Modified TINGZ
%       -add commetns

objimg=tz_objimg(obj,[]);

bwobjimg=objimg>0;

Theta = aw_angle(bwobjimg)*180./pi;
    
Rot = imrotate(objimg, Theta, 'bilinear');

[imgaxis,axln,dists,borders]=tz_imaxis(Rot);
[imgaxis2,axln2,dists2,borders2]=tz_imaxis(Rot');

if size(axln2,1)>size(axln,1)
    objimg=Rot';

else
    objimg=Rot;
end

%bwRot=Rot>0;
[ox,oy]=find(objimg>0);
cy=mean(oy);
if cy<size(objimg,1)/2
    objimg=flipud(objimg);
end

projy=sum(objimg,2);
holds=find(projy>0);
projx=sum(objimg,1);
holdsx=find(projx>0);

objimg=objimg(min(holds):max(holds),:);
objimg=objimg(:,min(holdsx):max(holdsx));

if isempty(imgsize)
    imgsize=size(objimg);
else
    objimg=imresize(objimg,[imgsize(1),size(objimg,2)],'bilinear');
end

img=zeros(imgsize);

objimg=[zeros(1,size(objimg,2));objimg;zeros(1,size(objimg,2))];

for i=1:imgsize(1)
    cs=objimg(i+1,:);
    holds=find(cs>0);
    if ~isempty(holds)
        if size(holds)==1
            img(i,:)=cs(holds);
        else
            subimg=objimg(i:i+2,min(holds):max(holds));
            imshow(subimg,[]);
            drawnow
            subimg=imresize(subimg,[3,imgsize(2)],'bilinear');
            img(i,:)=subimg(2,:);
        end
    end
end

% for i=1:imgsize(1)
%     cs=objimg(i,:);
%     holds=find(cs>0);
%     
%     if size(holds)==1
%         img(i,:)=objimg(i,holds);
%     else
%         ns=(0:imgsize(2)-1)/(imgsize(2)-1)*(max(holds)-min(holds))+min(holds);
%         img(i,:)=interp1(min(holds):max(holds),objimg(i,min(holds):max(holds)),ns,'spline');
%     end
% end
