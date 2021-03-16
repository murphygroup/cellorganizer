function [Circle_Center_X,Circle_Center_Y,R,min_tmp]=GetCircle_ByCenter_Distance(circleParaXYR,Circle_Center_X,Circle_Center_Y)
for k=1:size(circleParaXYR,1)
      tmp(k)=( (circleParaXYR(k,2)-Circle_Center_X)^2 + (circleParaXYR(k,1)-Circle_Center_Y)^2 );
end
[min_tmp,idex]=min(sqrt(tmp));
Circle_Center_X=circleParaXYR(idex,2);
Circle_Center_Y=circleParaXYR(idex,1);
R=circleParaXYR(idex,3);