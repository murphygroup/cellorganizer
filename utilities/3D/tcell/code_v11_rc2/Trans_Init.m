function [Point]=Trans_Init( DIC,GFP,lunkuo_of_ImgInput1,InitPoint,ZhiXin)
    for i=1:length(lunkuo_of_ImgInput1)
        %求出两点连线的直线
        k=(InitPoint(4)-InitPoint(2))/(InitPoint(3)-InitPoint(1));
        b=InitPoint(2)-k*InitPoint(1);
        TempLeft=[InitPoint(1),InitPoint(2)];
        TempRight=[InitPoint(3),InitPoint(4)];       
        center=[(InitPoint(1)+InitPoint(3))/2,(InitPoint(4)+InitPoint(2))/2];
        
        vector=ZhiXin(i,:)-center;
       Dis_Center=sqrt((vector(1)^2+vector(2)^2));
        vectorper=vector/Dis_Center;
        %求出轮廓与直线相交的点得到深度
        temp_k=(center(2)-ZhiXin(i,1))/(center(1)-ZhiXin(i,2));
        temp_b=center(2)-temp_k*center(1);
        y=lunkuo_of_ImgInput1{i,1}(:,1)*temp_k+temp_b*ones(size(lunkuo_of_ImgInput1{i,1},1),1);
        Piancha=y-lunkuo_of_ImgInput1{i,1}(:,2);
        [num,mm]=sort(Piancha);
        %找出5个点离两连线中间点最近就为轮廓顶点
        All_dis2=[];
         for j=1:5
            dis2=sqrt((ZhiXin(i,1)-lunkuo_of_ImgInput1{i,1}(mm(j),1))^2+(ZhiXin(i,2)-lunkuo_of_ImgInput1{i,1}(mm(j),2))^2);
            All_dis2=[All_dis2,dis2];
         end
        [num,m]=min(All_dis2);
        LuoKuo_Point=[lunkuo_of_ImgInput1{i,1}(mm(m),1),lunkuo_of_ImgInput1{i,1}(mm(m),2)];

        Temp_Gfp=zeros(size(GFP{i,1}));
        for j=1:size(lunkuo_of_ImgInput1{i,1},1)
            Temp_Gfp(floor(lunkuo_of_ImgInput1{i,1}(j,1)),floor(lunkuo_of_ImgInput1{i,1}(j,2)))=GFP{i,1}(floor(lunkuo_of_ImgInput1{i,1}(j,1)),floor(lunkuo_of_ImgInput1{i,1}(j,2)));
        end
        %如果轮廓上的点到质心的距离大于直线中点到质心的距离就沿质心和轮廓上的点（靠近大细胞的点）连线平移
        Dis_LuoKuo=sqrt((LuoKuo_Point(1)-center(1))^2+(LuoKuo_Point(2)-center(2))^2);   
        %如果轮廓上的点到质心小于直线中点到质心的距离就沿中点和质心连线向小细胞方向平移
         if  Dis_LuoKuo < Dis_Center      
            Shendu=Dis_LuoKuo-1/4*(Dis_LuoKuo);
            All_Avg_Gray=[];       
            while Dis_Center>= Shendu
                    TempLeft=TempLeft+vectorper;
                    TempRight=TempRight+vectorper;
%                     Temp_Gfp=zeros(GFP{i,1});
%                     for j=1:size(lunkuo_of_ImgInput1{i,1},1)
%                         Temp_Gfp(lunkuo_of_ImgInput1{i,1}(j,1),lunkuo_of_ImgInput1{i,1}(j,2))=GFP{i,1}(lunkuo_of_ImgInput1{i,1}(i,1),lunkuo_of_ImgInput1{i,1}(i,2));
%                     end
                    TempLine_k=(TempLeft(2)-TempRight(2))/(TempLeft(1)-TempRight(1));
                    TempLine_b=TempRight(2)-TempRight(1)*TempLine_k;
                    Avg_Gray=BigGray( Temp_Gfp,TempLine_k,TempLine_b,TempLeft,TempRight );
                    center=(TempLeft+TempRight)/2;
                    vector=ZhiXin(i,:)-center;
                    Dis_Center=sqrt((vector(1)^2+vector(2)^2));
                    All_Avg_Gray=[All_Avg_Gray,Avg_Gray];
            end
         else
                vector=LuoKuo_Point-center; 
                vectorper=vector/(sqrt(vector(1)^2+vector(2)^2));
                Shendu=Dis_LuoKuo-1/4*(Dis_LuoKuo);
                All_Avg_Gray=[];
                while Dis_Center < Shendu
                    TempLeft=TempLeft+vectorper;
                    TempRight=TempRight+vectorper;

                    TempLine_k=(TempLeft(2)-TempRight(2))/(TempLeft(1)-TempRight(1));
                    TempLine_b=TempRight(2)-TempRight(1)*TempLine_k;
                    Avg_Gray=BigGray( Temp_Gfp,TempLine_k,TempLine_b,TempLeft,TempRight );
                    center=(TempLeft+TempRight)/2;
                    vector=LuoKuo_Point-center;
                    Dis_Center=sqrt((vector(1)^2+vector(2)^2));
                    All_Avg_Gray=[All_Avg_Gray,Avg_Gray];
                end
         end
         [num,mmm]=max(All_Avg_Gray);
         TempLeft=TempLeft-(length(All_Avg_Gray)-mmm+1)*vectorper;
         TempRight=TempRight-(length(Avg_Gray)-mmm+1)*vectorper;
         Point{i,1}=[TempLeft(1),TempLeft(2),TempRight(1),TempRight(2)];
    end
end

