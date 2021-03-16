function [ZhixinX,ZhixinY,lunkuo_of_ImgInput1,Img]=InitialSmallCell(img,width,Excel_cp)
% Select image
img=medfilt2(img);
% img2=imcrop(img,[Excel_cp(1)-width/2,Excel_cp(2)-width/2,width,width]);
[img2,Position]=Imcrop(img,Excel_cp(1),Excel_cp(2),width);
% Otsu
level=graythresh(img2);BW=im2bw(img2,level);
BW=medfilt2(BW);
se = strel('disk',1);
BW = imclose(BW,se);
% Fill in the inner hole area
Img=imfill(BW,'holes');
% Find the connected domain and label it
[L,num] = bwlabel(Img,8);
% Calculate the centroid of each connected domain
distance=10000;
for j=1:num
    [m,n]=find(L==j);
    ZhiXinX=sum(n)/length(n);
    ZhiXinY=sum(m)/length(m);     
    %Calculate the distance from the center of mass to the upper left corner
    dis=sqrt((width/2-ZhiXinX)^2+(width/2-ZhiXinY)^2);
    distance=[distance;dis];
end
% Calculate the distance from the center of mass to the upper left corner, which is the shortest connected domain
distance=distance(2:end);
[min_num,ix]=min(distance);
% Find the point at which all values are the boundary, and the rest is 0 (black background)
[f_m,f_n]=find(L==ix);
[Im,In]=size(Img);
for j=1:Im
    for j_n=1:In
        Img(j,j_n)=0;
    end
end
for j=1:length(f_m)
    Img(f_m(j),f_n(j))=1;
end
% edge detection
edge_BW=edge(Img,'sobel');
[m,n]=find(edge_BW==1); 
 m=m+Position(2);n=n+Position(1);
lunkuo_of_ImgInput1=[n,m];     % This contour coordinate is in the small figure that draws down below, begin to restore to inside big graph next

% The centroid coordinates of the target region are obtained
ZhixinY=sum(m)/length(m);
ZhixinX=sum(n)/length(n);

