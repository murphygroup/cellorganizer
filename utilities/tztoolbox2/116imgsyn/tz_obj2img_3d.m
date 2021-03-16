function img=tz_obj2img_3d(objects,imgsize)
%TZ_OBJ2IMG_3D Obsolete.

%function img=tz_obj2img_3d(objects)
%
%OVERVIEW:
%   conver 3d objects to 3d images
%PARAMETERS:
%   objects - objects, cell array
%   imgsize - 3d image size
%RETURN:
%   img - 3d matrix
%DESCRIPTION:
%
%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   01-NOV-2004 Modified TINGZ
%       - add comments
%       - deal with multiple objects
%

Warning(tz_genmsg('of','tz_obj2img_3d','tz_objs2img'));

img=zeros(imgsize);

for i=1:length(objects)
    idx=tz_index2order(objects{i}(:,1:3),imgsize);
    idx(idx==0)=[];
    intensity=objects{i}(idx);
    img(idx)=intensity;
end
