function [Circle_LunKuo,Center,CircleR,Circle_Center_X,Circle_Center_Y]=InitialBigCell2(Excel_lp,Excel_rp,r,ZhixinX,ZhixinY)
syms a b
f1=(Excel_lp(1)-a)^2+(Excel_lp(2)-b)^2-r^2;
f2=(Excel_rp(1)-a)^2+(Excel_rp(2)-b)^2-r^2;
[a,b]=vpasolve(f1,f2,a,b);
A=double([a,b]);
result1=(A(1,1)-ZhixinX)^2+(A(1,2)-ZhixinY)^2;
result2=(A(2,1)-ZhixinX)^2+(A(2,2)-ZhixinY)^2;
if result1>result2
    result=A(1,:);
else
    result=A(2,:);
end
result=round(result);
Center=[result(1),result(2)];
CircleR=r;
theta = 0:0.05:6.28;
xx = result(1) + r*cos(theta);
yy = result(2) + r*sin(theta);
Circle_LunKuo=[xx',yy'];
Circle_Center_Y=result(2) ;
Circle_Center_X=result(1);




