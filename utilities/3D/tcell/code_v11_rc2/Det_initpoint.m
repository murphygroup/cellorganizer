function [ JD ] = Det_initpoint( initpoint, DIC,Superposition,Finally)
Excel_lp=[initpoint(1),initpoint(2)];
Excel_rp=[initpoint(3),initpoint(4)];
r=20;Start=initpoint(end-1);
Excel_cp=[(Excel_lp(1)+Excel_rp(1))/2,(Excel_lp(2)+Excel_rp(2))/2];
width1=34;width2=60;cnt=1;
[ZhixinX,ZhixinY,lunkuo_of_ImgInput1,Img]=InitialSmallCell(Superposition{Start,1},width1,Excel_cp);
target_lunkuo=lunkuo_of_ImgInput1;
lunkuo_of_ImgInput{Start,1}=lunkuo_of_ImgInput1;
ZhixinX1{Start,1}=ZhixinX;
ZhixinY1{Start,1}=ZhixinY;
 JiaoDian_Point{Start,1}=[Excel_lp;Excel_rp];
[Circle_LunKuo,Center,CircleR,Circle_Center_X,Circle_Center_Y]=InitialBigCell2(Excel_lp,Excel_rp,r,ZhixinX,ZhixinY);
 Need=[JiaoDian_Point{Start,1}(1,1),JiaoDian_Point{Start,1}(1,2),JiaoDian_Point{Start,1}(2,1),JiaoDian_Point{Start,1}(2,2),Start];
% figure(1)
% imshow(DIC{Start,1});
% hold on
% plot(JiaoDian_Point{Start,1}(1,1),JiaoDian_Point{Start,1}(1,2),'r*',JiaoDian_Point{Start,1}(2,1),JiaoDian_Point{Start,1}(2,2),'r*');

 for j = Start+1:Finally
        %Small cell segmentation
        [ZhixinX,ZhixinY,lunkuo_of_ImgInput1,Img]=SmallCellTracking(Superposition{j,1},ZhixinX,ZhixinY,width1,Img);
        if (abs(ZhixinX-ZhixinX1{j-1,1})>15) &&(abs(ZhixinY-ZhixinY1{j-1,1})>15)
            lunkuo_of_ImgInput1=lunkuo_of_ImgInput{j-1,1};
            ZhixinX=ZhixinX1{j-1,1};
            ZhixinY=ZhixinY1{j-1,1};
        end
        lunkuo_of_ImgInput{j,1}=lunkuo_of_ImgInput1;
        ZhixinX1{j,1}=ZhixinX;
        ZhixinY1{j,1}=ZhixinY;
        %Great cell segmentation
        [Center,CircleR,Circle_LunKuo,Circle_Center_X,Circle_Center_Y]=BigCellTracking(DIC{j,1},Center,CircleR,Circle_Center_X,Circle_Center_Y,width2);  
        
%         image1 = imread([file_path3,strcat(int2str(j),'.png')],'png');
        for i_q=j-1:-1:Start
            if ~isempty(JiaoDian_Point{i_q,1})
                JiaoDian=JiaoDian_Point{i_q,1};
                break;
            end
        end
        if isempty(lunkuo_of_ImgInput1)
            for i_q=j-1:-1:Start
                if ~isempty(lunkuo_of_ImgInput{i_q,1})
                    lunkuo_of_ImgInput1=lunkuo_of_ImgInput{i_q,1};
                    ZhixinX=ZhixinX1{i_q,1};
                    ZhixinY=ZhixinY1{i_q,1};
                    break;
                end
            end
        end
        %Calculate intersection point    
        [JiaoDian_Point{j,1}]=DrawPointsInLoop(j,lunkuo_of_ImgInput1,Circle_LunKuo,Center,CircleR,JiaoDian, Superposition{j, 1});
%         figure(1)
%     imshow(DIC{j,1});hold on
%     plot(lunkuo_of_ImgInput1(:,1),lunkuo_of_ImgInput1(:,2),'r.');hold on
%     plot(Circle_LunKuo(:,1),Circle_LunKuo(:,2),'g.');
%     
%     if ~isempty(JiaoDian_Point{j,1})
%             plot(JiaoDian_Point{j,1}(1,1),JiaoDian_Point{j,1}(1,2),'r*',JiaoDian_Point{j,1}(2,1),JiaoDian_Point{j,1}(2,2),'r*');
%     end
        %Determine the right and left intersection points
        if ~isempty(JiaoDian_Point{j,1})
            Need=[Need;JiaoDian_Point{j,1}(1,1),JiaoDian_Point{j,1}(1,2),JiaoDian_Point{j,1}(2,1),JiaoDian_Point{j,1}(2,2),j];
            JD=[JiaoDian_Point{j,1}(1,1),JiaoDian_Point{j,1}(1,2),JiaoDian_Point{j,1}(2,1),JiaoDian_Point{j,1}(2,2),j];
        else
            Need=[Need;0,0,0,0,j];
            JD=[0,0,0,0,j];
        end 
end
end

