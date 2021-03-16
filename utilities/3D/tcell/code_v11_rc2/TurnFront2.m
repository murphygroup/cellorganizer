function [Need]=TurnFront2(Overlay,DIC,GFP, Start,Excel_lp, Excel_rp,inferred_annotation_directory,param)

% SaveName = strcat(inferred_annotation_directory,'Result222.csv');%Excel path and name to save the final result
% Processdir=strcat(inferred_annotation_directory,'/',num2str(k));
% if ~exist(Processdir, 'dir')
%     mkdir(Processdir);
% end
Processdir=inferred_annotation_directory;
JiaoDian_Point{Start,1}=[Excel_lp;Excel_rp];
r=20;%Empirical radius of great circle
Excel_cp=[(Excel_lp(1)+Excel_rp(1))/2,(Excel_lp(2)+Excel_rp(2))/2];
width1=34;width2=60;cnt=1;
% Small cell segmentation initialization
[ZhixinX,ZhixinY,lunkuo_of_ImgInput1,Img]=InitialSmallCell(GFP{Start,1},width1,Excel_cp);
target_lunkuo=lunkuo_of_ImgInput1;
lunkuo_of_ImgInput{Start,1}=lunkuo_of_ImgInput1;
ZhixinX1{Start,1}=ZhixinX;
ZhixinY1{Start,1}=ZhixinY;
% Large cell segmentation initialization
[Circle_LunKuo,Center,CircleR,Circle_Center_X,Circle_Center_Y]=InitialBigCell2(Excel_lp,Excel_rp,r,ZhixinX,ZhixinY);
% Determine the first frame of small cells in the location of large cells
% and distinguish the left and right point
if isempty(JiaoDian_Point{Start,1})
    Need=[0,0,0,0,Start];
else
    Need=[JiaoDian_Point{Start,1}(1,1),JiaoDian_Point{Start,1}(1,2),JiaoDian_Point{Start,1}(2,1),JiaoDian_Point{Start,1}(2,2),Start];%��ʼ֡excel���?
end

p_all=[];ZhiXin_all=[];xmin_all=[];xmax_all=[];ymin_all=[];ymax_all=[];All_Circle_Center=[];
for j = Start-1:-1:max(Start-2,1)
    %Small cell segmentation
    [ZhixinX,ZhixinY,lunkuo_of_ImgInput1,Img, Position]=SmallCellTracking(GFP{j,1},ZhixinX,ZhixinY,width1,Img);
    Superposition{abs(j-Start),1}=GFP{j,1};
    
    if abs(ZhixinX-ZhixinX1{j+1,1})>15 ||abs(ZhixinY-ZhixinY1{j+1,1})>15
        lunkuo_of_ImgInput1=lunkuo_of_ImgInput{j+1,1};
        ZhixinX=ZhixinX1{j+1,1};
        ZhixinY=ZhixinY1{j+1,1};
    end
    lunkuo_of_ImgInput{j,1}=lunkuo_of_ImgInput1;
    ZhixinX1{j,1}=ZhixinX;
    ZhixinY1{j,1}=ZhixinY;
    %Great cell segmentation
    [Center,CircleR,Circle_LunKuo,Circle_Center_X,Circle_Center_Y]=BigCellTracking2(DIC{j,1},j,Center,CircleR,Circle_Center_X,Circle_Center_Y,width2,GFP{j,1},Processdir, param);
    
    image1 = DIC{j,1};
    
    GFP_raw = GFP{j,1};
    GFP_raw_crop = GFP_raw(round(Position(2)) + (1:size(Img, 1))-1, round(Position(1)) + (1:size(Img, 2))-1);
    GFP_raw_seg = double(GFP_raw_crop) .* double(Img);
    stats = regionprops(Img);
    centroid = stats.Centroid;
    left_right_end_points = [Excel_lp, Excel_rp] - [Position(1:2), Position(1:2)];
    centroid_apc = [Circle_Center_X, Circle_Center_Y] - Position(1:2);
    Point = Trans_Init_1( GFP_raw_seg, Img, centroid, centroid_apc, left_right_end_points);
    if ~isempty(Point)
        Point_1 = Point + [Position(1:2), Position(1:2)];
    else
        Point_1 = [0, 0, 0, 0];
    end
    if false
        figure, imshow(GFP_raw_crop, []);
    end
        
    Need = [Need; Point_1, j];

end
if size(Need,1)==3
    InitPoint=Need(1,:);
    Need=Need(2:3,:);
    if(Need(1,1)==0) && (Need(2,1)==0)
        Need=[];
    end
elseif  size(Need,1)==2
      Need=Need(2,:);
else
     Need=[];
end
