function [ZhixinX,ZhixinY,lunkuo_of_ImgInput1,Img, Position]=SmallCellTracking(ImgInput2,ZhixinX,ZhixinY,width,Img)
% temp=imcrop(ImgInput2,[ZhixinX-width/2,ZhixinY-width/2,width,width]);
[temp,Position]=Imcrop(ImgInput2,ZhixinX,ZhixinY,width);
level=graythresh(temp);BW=im2bw(temp,level);

se = strel('disk',1);
BW = imclose(BW,se);
bb=imfill(BW,'holes');
[L,num] = bwlabel(bb,8);
%small area of connected will be removed
Region_Area=regionprops(L,'Area');
if (~isempty(Region_Area))
    for i=1:num
        if(Region_Area(i).Area<30)            
            [m,n]=find(L==i);
            bb(m,n)=0;
        end
    end
end
[L,num] = bwlabel(bb,8);
if num==0
    lunkuo_of_ImgInput1=[];
    Relative=1;
    return;
end
Zhi_Xin_L=regionprops(L,'Centroid');
ttmp=L;
theta=0:30:360;
area_of_every=zeros(num,length(theta));
all_max_area=zeros(1,length(theta));
all_max_location=zeros(1,length(theta));
for kk=0:30:360
    B = imrotate(Img,kk,'nearest','crop');
    for i=1:num
        L(ttmp~=i)=0;
        L(ttmp==i)=255;
        L_moved=move_img_to_center(L);
        L_moved(L_moved==-Inf)=0;
        area_of_every(i,kk/30+1)=sum(sum(B.*L_moved==255));
    end  
     [all_max_area(kk/30+1),all_max_location(kk/30+1)]=max(area_of_every(:,kk/30+1));
end
    [max_area,max_location]=max(all_max_area);
    location=all_max_location(max_location);
    
    zhixin_of_min_area=Zhi_Xin_L(location).Centroid;
    %The centroid of the cell position, which is obtained by minimizing the area overlap, serves as a new centroid
    ZhixinX=zhixin_of_min_area(1)+Position(1);
    ZhixinY=zhixin_of_min_area(2)+Position(2);
    % A new centroid is used to extract small cells in the original drawing
%     img2=imcrop(ImgInput2,[ZhixinX-width/2,ZhixinY-width/2,width,width]);
    [img2,Position]=Imcrop(ImgInput2,ZhixinX,ZhixinY,width);
    if Position(1) > 500 || Position(2) > 444
        flag = 1;
    end
%     if param.verbose
%         imwrite(img2,strcat(Processdir_small,'/','02NewCrop.png'),'png');
%     end
    % Otsu
    level2=graythresh(img2);BW2=im2bw(img2,level2);
%     subplot(122),imshow(BW2)
%     figure(101),imshow(img44);
    % Fill in the inner hole area
    Img2=imfill(BW2,'holes');
    L2= bwlabel(Img2,8);
    stats = regionprops(L2);
    Ar = cat(1, stats.Area);
    ind = find(Ar ==max(Ar));%Find the label of the largest connected region
    if length(ind) > 1
		Centroids = cat(1, stats.Centroid);
        Maxima_centroids = Centroids(ind, :);
        [~, closest_ind] = min(dist2(Maxima_centroids, ([size(L2, 2), size(L2, 1)] + 1) ./ 2));
        ind = ind(closest_ind);
    end    
    Img2(L2~=ind)=0;%Set the rest of the area to 0
    % edge detection    
    Img=Img2;
    % imwrite(Img,strcat(Processdir_small,'\','03Newcell.png'),'png');
    edge_BW2=edge(Img2,'sobel');
%     if param.verbose
%         imwrite(Img,strcat(Processdir_small,'/','03Newcell.png'),'png');
%         imwrite(edge_BW2,strcat(Processdir_small,'/','04edge.png'),'png');
%     end
    [m2,n2]=find(edge_BW2==1);
    m2=m2+Position(2);n2=n2+Position(1);
    lunkuo_of_ImgInput1=[n2,m2];     % This contour coordinate is in the small figure that draws down below, begin to restore to inside big graph next
    ZhixinY=sum(m2)/length(m2);   
    ZhixinX=sum(n2)/length(n2);
    % If the upper cell and the cell transformation area of this frame does not exceed 200, it is considered that it does not coincide