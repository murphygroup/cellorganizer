function tz_generatecell(nucbody)
%TZ_GENERATECELL Unknown.

%function tz_generatecell(nucbody)

nucimg=tz_obj2img(nucbody,imgsize,{'2d','bn'});
nuc=bwperim(nucimg);

theta=tz_bwmajorangle(nucimg)*180/pi;
len=sqrt(sum(imgsize.^2));
celldist=[];
nucdist=[];
angles=0:da:360-da;
for a=angles
    ext=sin(a*pi/180)*50+100;
    pts=tz_getlinept2(nuccenter(1:2),a+theta,len);
    %pts2=tz_getlinept2(nuccenter(1:2),a+theta,len+ext);
    ps=tz_imgptspixel(nucedge,pts);
    intc=find(ps>0);
    if ~isempty(intc)
        pts2=tz_getlinept2(nuccenter(1:2),a+theta,intc(1)+ext);
        nucimg(pts(:,1),pts(:,2))=1;
    end
end

imshow(nucimg,[]);

