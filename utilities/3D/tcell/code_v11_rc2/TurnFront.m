function [ Need ] = TurnFront( Overlay,DIC,GFP, Start,Excel_lp, Excel_rp)
JiaoDian_Point{Start,1}=[Excel_lp;Excel_rp];
r=20;%Empirical radius of great circle
Excel_cp=[(Excel_lp(1)+Excel_rp(1))/2,(Excel_lp(2)+Excel_rp(2))/2];
width1=34;width2=60;
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
    Need=[JiaoDian_Point{Start,1}(1,1),JiaoDian_Point{Start,1}(1,2),JiaoDian_Point{Start,1}(2,1),JiaoDian_Point{Start,1}(2,2),Start];%��ʼ֡excel���
end
% image1 = Overlay{Start,1};
% figure(1);imshow(image1);hold on
% plot(target_lunkuo(:,1),target_lunkuo(:,2),'r.');
% hold on
% plot(Circle_LunKuo(:,1),Circle_LunKuo(:,2),'g.');
% hold on
% plot(JiaoDian_Point{Start,1}(1,1),JiaoDian_Point{Start,1}(1,2),'r*',JiaoDian_Point{Start,1}(2,1),JiaoDian_Point{Start,1}(2,2),'y*');
% hold on
% plot(ZhixinX,ZhixinY,'yo');

for j = Start-1:-1:1
    %Small cell segmentation
    [ZhixinX,ZhixinY,lunkuo_of_ImgInput1,Img]=SmallCellTracking(GFP{j,1},ZhixinX,ZhixinY,width1,Img);
    if abs(ZhixinX-ZhixinX1{j+1,1})>15 || abs(ZhixinY-ZhixinY1{j+1,1})>15
        lunkuo_of_ImgInput1=lunkuo_of_ImgInput{j+1,1};
        ZhixinX=ZhixinX1{j+1,1};
        ZhixinY=ZhixinY1{j+1,1};
    end
    lunkuo_of_ImgInput{j,1}=lunkuo_of_ImgInput1;
    ZhixinX1{j,1}=ZhixinX;
    ZhixinY1{j,1}=ZhixinY;
    %Great cell segmentation
    [Center,CircleR,Circle_LunKuo,Circle_Center_X,Circle_Center_Y]=BigCellTracking(DIC{j,1},Center,CircleR,Circle_Center_X,Circle_Center_Y,width2);
    
    image1 = Overlay{j,1};
%     figure(1);clf,imshow(image1);hold on
%     plot(Circle_Center_X,Circle_Center_Y,'ro');
%     hold on
    
    for i_q=j+1:1:Start
        if ~isempty(JiaoDian_Point{i_q,1})
            JiaoDian=JiaoDian_Point{i_q,1};
            break;
        end
    end
    if isempty(lunkuo_of_ImgInput1)
        for i_q=j+1:1:Start
            if ~isempty(lunkuo_of_ImgInput{i_q,1})
                lunkuo_of_ImgInput1=lunkuo_of_ImgInput{i_q,1};
                ZhixinX=ZhixinX1{i_q,1};
                ZhixinY=ZhixinY1{i_q,1};
                break;
            end
        end
    end
    %Calculate intersection point
    [JiaoDian_Point{j,1}]=DrawPointsInLoop(j,lunkuo_of_ImgInput1,Circle_LunKuo,Center,CircleR,JiaoDian,GFP{j,1});
    %Determine the right and left intersection points
    if ~isempty(JiaoDian_Point{j,1})
        Need=[Need;JiaoDian_Point{j,1}(1,1),JiaoDian_Point{j,1}(1,2),JiaoDian_Point{j,1}(2,1),JiaoDian_Point{j,1}(2,2),j];
%         plot(JiaoDian_Point{j,1}(1,1),JiaoDian_Point{j,1}(1,2),'r*',JiaoDian_Point{j,1}(2,1),JiaoDian_Point{j,1}(2,2),'y*');
        %         saveas(gcf,[directory,int2str(j)],'png');
    else
        Need=[Need;0,0,0,0,j];
    end
end
end

