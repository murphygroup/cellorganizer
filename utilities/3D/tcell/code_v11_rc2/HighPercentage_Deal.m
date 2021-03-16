function  JiaoDian_Point  = HighPercentage_Deal(Bx,By,JiaoDian_Point,p,Superposition)  

CenterX=(JiaoDian_Point(1,1)+JiaoDian_Point(2,1))/2;
CenterY=(JiaoDian_Point(1,2)+JiaoDian_Point(2,2))/2;
dis=sqrt((Bx-CenterX)^2+(By-CenterY)^2);
Vect=[(Bx-CenterX)/dis,(By-CenterY)/dis];
TempLeft=JiaoDian_Point(1,:);
TempRight=JiaoDian_Point(2,:);
TempLine_k=(TempLeft(2)-TempRight(2))/(TempLeft(1)-TempRight(1));
TempLine_b=TempRight(2)-TempRight(1)*TempLine_k;
Avg_Gray=BigGray( Superposition,TempLine_k,TempLine_b,TempLeft,TempRight );
while(1)
    TempLeft=TempLeft+Vect;
    TempRight=TempRight+Vect;
    TempLine_k=(TempLeft(2)-TempRight(2))/(TempLeft(1)-TempRight(1));
    TempLine_b=TempRight(2)-TempRight(1)*TempLine_k;
     syms x y 
    [x_solve,y_solve]=vpasolve(p(1)*x^2+p(2)*x*y+p(3)*y^2+p(4)*x+p(5)*y+p(6)==0,TempLine_k*x+TempLine_b==y,x,y);
    solutions2=double([x_solve,y_solve]);
    solutions2 = solutions2(imag(solutions2(:, 1)) == 0, :);
    
    if (size(solutions2) >= 2)
        dis2=sqrt((solutions2(1,1)-solutions2(2,1))^2+(solutions2(1,2)-solutions2(2,2))^2);
        if(dis2<10)
            break;
        else
           Avg_Gray=[Avg_Gray;BigGray( Superposition,TempLine_k,TempLine_b,solutions2(1,:),solutions2(2,:))];
        end   
    else
        return;
    end
end
[num,m]=max(Avg_Gray);
TempLeft=TempLeft-(length(Avg_Gray)-m+1)*Vect;
TempRight=TempRight-(length(Avg_Gray)-m+1)*Vect;

TempLine_k=(TempLeft(2)-TempRight(2))/(TempLeft(1)-TempRight(1));
TempLine_b=TempRight(2)-TempRight(1)*TempLine_k;
syms x y 
% [x_solve,y_solve]=solve(p(1)*x^2+p(2)*x*y+p(3)*y^2+p(4)*x+p(5)*y+p(6)==0,TempLine_k*x+TempLine_b==y,'Real',true);
 [x_solve,y_solve]=vpasolve(p(1)*x^2+p(2)*x*y+p(3)*y^2+p(4)*x+p(5)*y+p(6)==0,TempLine_k*x+TempLine_b==y,x,y,[-Inf Inf]);
solutions2=double([x_solve,y_solve]);
 solutions2 = solutions2(imag(solutions2(:, 1)) == 0, :);
if (size(solutions2) >= 2)
    JiaoDian_Point=solutions2;
else
    return;
end

