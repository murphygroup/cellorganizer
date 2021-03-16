function Percentage = Percentage_overlap( p,Sxmin,Sxmax,Symin,Symax,Bx,By,BR )
overlap_t=0;
t=0;
n=50000;
for i=1:n
  y1=round((Symax-Symin)*rand()+Symin);
  x1=round((Sxmax-Sxmin)*rand()+Sxmin);
  %The area of the overlapping part
  if  p(1)*x1^2+p(2)*x1*y1+p(3)*y1^2+p(4)*x1+p(5)*y1+p(6)>0 && (x1-Bx)^2+(y1-By)^2<(BR)^2
      overlap_t=overlap_t+1;
  end
  %The area of the ellipse
  if  p(1)*x1^2+p(2)*x1*y1+p(3)*y1^2+p(4)*x1+p(5)*y1+p(6)>0 
      t=t+1;
  end
end
Percentage=overlap_t/t;
end

