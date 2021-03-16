function Img=move_img_to_center(img)
[height,width]=size(img);
[xx,yy]=find(img~=0);
center_xx=sum(xx)/length(xx);
center_yy=sum(yy)/length(yy);
se=translate(strel(1),[fix( (height/2-center_xx) )   fix((width/2-center_yy))]);
Img=imdilate(img,se);