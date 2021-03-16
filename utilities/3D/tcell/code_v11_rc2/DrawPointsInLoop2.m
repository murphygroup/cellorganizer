function [JiaoDian_Point,cnt]=DrawPointsInLoop2(j,target_lunkuo,Circle_LunKuo,Center,CircleR,cnt,JiaoDian_Point_shang, Superposition, param)
% image1 = imread([file_path1,strcat(int2str(j),'.png')],'png');
% imshow(image1);hold on;
%  hold on; plot(target_lunkuo(:,1),target_lunkuo(:,2),'r.');% Draw a small cell outline
% Here's an elliptical fitting of small cells
[a,b,q,p,p2,X_center,Y_center,xmin,xmax,ymin,ymax,F]=Ellipse_Fitting3([target_lunkuo(:,1),target_lunkuo(:,2)], param);

if param.verbose
    hold on
    plot(Circle_LunKuo(:,1),Circle_LunKuo(:,2),'g.');% Draw a large cell outline
    hold on
end

syms x y 
x0=Center(1);y0=Center(2);R=CircleR;
% Compute the intersection of an ellipse and a circle
% [x_solve,y_solve]=solve(p(1)*x^2+p(2)*x*y+p(3)*y^2+p(4)*x+p(5)*y+p(6)==0, (x-x0)^2+(y-y0)^2== R^2);
[x_solve,y_solve]=vpasolve(p(1)*x^2+p(2)*x*y+p(3)*y^2+p(4)*x+p(5)*y+p(6)==0, (x-x0)^2+(y-y0)^2== R^2,x,y,[-Inf Inf]);
% Solve=solve(-p(1)*x^2-p(2)*x*y-p(3)*y^2-p(4)*x-p(5)*y-p(6)==0, (x-x0)^2+(y-y0)^2== R^2,'Real',true);
solutions2=double([x_solve,y_solve]);
x_solve=solutions2(:,1);y_solve=solutions2(:,2);
idex=find(imag(x_solve)==0);                        
if ~isempty(idex)
    Real_x_solve=x_solve(idex);
    Real_y_solve=y_solve(idex);
%     hold on;
%     plot(Real_x_solve(1), Real_y_solve(1), 'r.', 'MarkerSize', 25);
%     hold on
%     plot(Real_x_solve(2), Real_y_solve(2), 'y.', 'MarkerSize', 25);
%     hold on
%     text(sum(Real_x_solve)/length(Real_x_solve)-20,sum(Real_y_solve)/length(Real_y_solve)-20,num2str(cnt),'FontSize',18);
    One_One=abs(Real_x_solve(1)-JiaoDian_Point_shang(1,1))+abs(Real_y_solve(1)-JiaoDian_Point_shang(1,2));
    One_Two=abs(Real_x_solve(1)-JiaoDian_Point_shang(2,1))+abs(Real_y_solve(1)-JiaoDian_Point_shang(2,2));
    Two_One=abs(Real_x_solve(2)-JiaoDian_Point_shang(1,1))+abs(Real_y_solve(2)-JiaoDian_Point_shang(1,2));
    Two_Two=abs(Real_x_solve(2)-JiaoDian_Point_shang(2,1))+abs(Real_y_solve(2)-JiaoDian_Point_shang(2,2));
    Sum_dis=One_One+Two_Two;
    Sum_dis2=One_Two+Two_One;
    if Sum_dis2<Sum_dis
        temp=Real_x_solve(1);
        Real_x_solve(1)=Real_x_solve(2);
        Real_x_solve(2)=temp;
        temp=Real_y_solve(1);
        Real_y_solve(1)=Real_y_solve(2);
        Real_y_solve(2)=temp;
    end
    cnt=cnt+1;
    JiaoDian_Point=[Real_x_solve,Real_y_solve];
    % judge the percentage of overlap
    Percentage = Percentage_overlap(p, xmin, xmax, ymin, ymax, x0, y0, R);
    if Percentage > 0.4
        JiaoDian_Point = HighPercentage_Deal(x0, y0, JiaoDian_Point, p, Superposition);
    end
else
    JiaoDian_Point=[];
end
% text(350,430,['StartFrame:' int2str(j)],'FontSize',15,'Color','BLACK');
