function [Need]=TurnFront1(Overlay,DIC,GFP, Start,Excel_lp, Excel_rp,inferred_annotation_directory,param)

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
% Processdir_big=strcat(Processdir,'/Big/',num2str(Start));
% if ~exist(Processdir_big, 'dir')
%     mkdir(Processdir_big);
% end
% image1 = DIC{Start,1};
% figure(1);imshow(image1);hold on
% plot(target_lunkuo(:,1),target_lunkuo(:,2),'r.');
% hold on
% plot(Circle_LunKuo(:,1),Circle_LunKuo(:,2),'g.');
% hold on
% plot(JiaoDian_Point{Start,1}(1,1),JiaoDian_Point{Start,1}(1,2),'r*',JiaoDian_Point{Start,1}(2,1),JiaoDian_Point{Start,1}(2,2),'y*');
% hold on
% plot(ZhixinX,ZhixinY,'yo');
% saveas(gcf,[strcat(Processdir_big,'\'),'Result.png'],'png');

p_all=[];ZhiXin_all=[];xmin_all=[];xmax_all=[];ymin_all=[];ymax_all=[];All_Circle_Center=[];
for j = Start-1:-1:max(Start-2,1)
    %Small cell segmentation
    [ZhixinX,ZhixinY,lunkuo_of_ImgInput1,Img]=SmallCellTracking(GFP{j,1},ZhixinX,ZhixinY,width1,Img);
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
     if param.verbose
             figure(1);clf,imshow(image1);hold on
             plot(Circle_Center_X,Circle_Center_Y,'ro');
             hold on
     end
    All_Circle_Center=[All_Circle_Center;Circle_Center_X,Circle_Center_Y];
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
    [JiaoDian_Point{j,1},p,ZhiXin,xmin,xmax,ymin,ymax]=DrawPointsInLoop3(j,lunkuo_of_ImgInput1,Circle_LunKuo,Center,CircleR,JiaoDian,GFP{j,1},param);
    p_all=[p_all;p'];ZhiXin_all=[ZhiXin_all;ZhiXin];xmin_all=[xmin_all;xmin];xmax_all=[xmax_all;xmax];ymin_all=[ymin_all;ymin];ymax_all=[ymax_all;ymax];
    LuoKuo{Start-j,1}=lunkuo_of_ImgInput1;
    GFP_Need{Start-j,1}=GFP{j,1};
    DIC_Need{Start-j,1}=DIC{j,1};
    %Determine the right and left intersection points
    if ~isempty(JiaoDian_Point{j,1})
        Need=[Need;JiaoDian_Point{j,1}(1,1),JiaoDian_Point{j,1}(1,2),JiaoDian_Point{j,1}(2,1),JiaoDian_Point{j,1}(2,2),j];
           if param.verbose
                    plot(JiaoDian_Point{j,1}(1,1),JiaoDian_Point{j,1}(1,2),'r*',JiaoDian_Point{j,1}(2,1),JiaoDian_Point{j,1}(2,2),'y*');
           end
        %         saveas(gcf,[directory,int2str(j)],'png');
    else
        Need=[Need;0,0,0,0,j];
    end
    Processdir_result=strcat(Processdir,'/result');
    if param.verbose
        if ~exist(Processdir_result, 'dir')
            mkdir(Processdir_result);
        end
        saveas(gcf,[strcat(Processdir_result,'/'), int2str(j), '.png'],'png');
    end
end
if size(Need,1)==3
    InitPoint=Need(1,:);
    Need=Need(2:3,:);
    %两个没交点
    if(Need(1,1)==0) && (Need(2,1)==0)
        
       Point= Trans_Init( DIC,GFP_Need,LuoKuo,InitPoint,ZhiXin_all);
       Need=[Point{1,1} Start-2;Point{2,1}, Start-1];
%          ResultPoint=Find_long_short( p_all,ZhiXin_all,xmin_all,xmax_all,All_Circle_Center,Superposition);
%          Frame=Need(:,end);
%          Need=[];
%          for i=1:length(ResultPoint)
%             Need=[Need;ResultPoint{i,1}(1,1),ResultPoint{i,1}(1,2),ResultPoint{i,1}(2,1),ResultPoint{i,1}(2,2),Frame(i)];
%          end
%          
% 		 if param.verbose
%          figure(1),imshow(DIC{Need(1,5),1});
%          hold on;plot(Need(1,1),Need(1,2),'y*');
%          hold on;plot(Need(1,3),Need(1,4),'y*');
%          if ~exist(Processdir_result, 'dir')
%             mkdir(Processdir_result);
%          end
%          saveas(gcf,[strcat(Processdir_result,'/'), int2str(Need(1,5)), '.png'],'png');
%          figure(1),imshow(DIC{Need(2,5),1});
%          hold on;plot(Need(2,1),Need(2,2),'y*');
%          hold on;plot(Need(2,3),Need(2,4),'y*');
%          saveas(gcf,[strcat(Processdir_result,'/'), int2str(Need(2,5)), '.png'],'png');
%          end
% %          %最前一帧有交点，后面一帧没交点
%      elseif  (Need(1,1)~=0) && (Need(2,1)==0)
            
    end
    Temp_need=Need(1,:);
    Need(1,:)=Need(2,:);
    Need(2,:)=Temp_need;
elseif  size(Need,1)==2
      Need=Need(2,:);
else
     Need=[];
end
% Need = Interpolation_3(Need); %
%Save result
% FormHead={'Lx','Ly','Rx','Ry','frame'};
% Need=num2cell(Need);
% Result = cell2table(Need,'VariableNames',FormHead);
% writetable(Result,SaveName,'WriteRowNames',true);

% All_Need=[FormHead;All_Need];
% xlswrite(SaveName,All_Need);
