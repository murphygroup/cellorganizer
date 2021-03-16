function  ResultPoint=Find_long_short( p,ZhiXin,xmin,xmax,Circle_Center,Superposition)
%To find the equations of the long axis and the short axis£¨four top points£©
%the nearest two points from the large cell will be moved in the direction of large cells
%Find the two top points of each ellipse's long axis and minor axis
for i=1:length(xmin)
    x=xmin(i):0.1:xmax(i);
    Point=[];
    for j=1:length(x)
        A=p(i,3);B=p(i,2)*x(j)+p(i,5);C=p(i,1)*x(j)^2+p(i,4)*x(j)+p(i,6);
        delt=B^2-4*A*C;
        if delt>0
            y1=(-B+sqrt(delt))/(2*A);
            y2=(-B-sqrt(delt))/(2*A);
        elseif delt==0
            y2=-B/(2*A);
            y1=y2;
        else
            y1=[];
            y2=[];
        end
        if ~isempty(y1)
            Point=[Point;x(j),y1;x(j),y2];
        else
            Point=Point;
        end
    end  
    dis=[];
    for j=1:size(Point,1)
        dis(j)=(ZhiXin(i,1)-Point(j,1))^2+(ZhiXin(i,2)-Point(j,2))^2;
    end
    [num1,short_cnt]=min(dis);
    [num2,big_cnt]=max(dis);
    short_Poistion=[Point(short_cnt,1),Point(short_cnt,2)];
    long_Poistion=[Point(big_cnt,1),Point(big_cnt,2)];
    %another point
    Other_short_Poistion=[ZhiXin(i,1)*2-short_Poistion(1),ZhiXin(i,2)*2-short_Poistion(2)];
    Other_long_Poistion=[ZhiXin(i,1)*2-long_Poistion(1),ZhiXin(i,2)*2-long_Poistion(2)];
    Top_Point=[short_Poistion;long_Poistion;Other_short_Poistion;Other_long_Poistion];
	%The two symmetric points that are not nearest to the center of mass of the great circle are intersection points, and then translation
    Dis=[];
    for j=1:size(Top_Point,1)
        Dis(j)=(Circle_Center(i,1)-Top_Point(j,1))^2+(Circle_Center(i,2)-Top_Point(j,2))^2;
    end
    [num,row]=min(Dis);
    Small_poistion=Top_Point(row,:);
    Other_Point=[ZhiXin(i,1)*2-Top_Point(row,1),ZhiXin(i,2)*2-Top_Point(row,2)];
    Top_Point(row,:)=[];
    for j=1:size(Top_Point,1)-1
        if ~isempty(find(Top_Point(j,:)==Other_Point))
            Top_Point(j,:)=[];
            break;
        end
    end
    JiaoDian_Point=Top_Point;
        ResultPoint{i,1}  = HighPercentage_Deal2(ZhiXin(i,1),ZhiXin(i,2),Small_poistion,JiaoDian_Point,p(i,:),Superposition{i,1});
end

