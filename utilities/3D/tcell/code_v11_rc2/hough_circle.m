function [hough_space,hough_circle,para] = hough_circle(BW,step_r,step_angle,r_min,r_max,p)  
  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% input  
% BW:two value image£»  
% step_r:Circle radius step 
% step_angle:angle step, in radians 
% r_min:Minimum radius of circle  
% r_max:maximum circle radius
% p:threshold, the number between 0 and 1. By adjusting this value, the center of the circle and the radius of the circle in the graph can be obtained 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% output  
% hough_space:The argument space, h(a, b, r), represents the number of points in the circle (a, b) with a radius of r
% hough_circl:Two valued image, detected circle 
% para:The center and radius of the circle to be detected 
  
circleParaXYR=[];  
para=[];  
  
[m,n] = size(BW);  
size_r = round((r_max-r_min)/step_r)+1; 
size_angle = round(2*pi/step_angle);  
  
hough_space = zeros(m,n,size_r);  
  
[rows,cols] = find(BW); 
ecount = size(rows);
  
% Hough transform
% corresponds the image space (x, y) to the parameter space (a, b, r)
  
% a = x-r*cos(angle)  
% b = y-r*sin(angle)  
  
for i=1:ecount  
    for r=1:size_r %Radius step number
        for k=1:size_angle %Divide the circle in a certain arc 
            a = round(rows(i)-(r_min+(r-1)*step_r)*cos(k*step_angle));  
            b = round(cols(i)-(r_min+(r-1)*step_r)*sin(k*step_angle));  
            if(a>0&a<=m&b>0&b<=n)  
            hough_space(a,b,r) = hough_space(a,b,r)+1;%h(a,b,r) h(centers and radii)  
            end  
        end  
    end  
end  
  
  
% Search for clustering points beyond the threshold. For the detection of multiple circles, the threshold should be set a little smaller! By adjusting this value, we can find the center and radius of all circles  
max_para = max(max(max(hough_space)));%The return value is the maximum value of this matrix  
index = find(hough_space>=max_para*p);%In a matrix, you want to find the location in which the number is greater than max_para*p 
length = size(index);%The number that conforms to the threshold 
hough_circle = false(m,n);  
%hough_circle = zeros(m,n);  
%Find the radius and center of the circle by position 
for i=1:ecount  
    for k=1:length  
        par3 = floor(index(k)/(m*n))+1;  
        par2 = floor((index(k)-(par3-1)*(m*n))/m)+1;  
        par1 = index(k)-(par3-1)*(m*n)-(par2-1)*m;  
        if((rows(i)-par1)^2+(cols(i)-par2)^2<(r_min+(par3-1)*step_r)^2+5&...  
                (rows(i)-par1)^2+(cols(i)-par2)^2>(r_min+(par3-1)*step_r)^2-5)  
              hough_circle(rows(i),cols(i)) = true;   %Detected circle 
        end  
    end  
end                 
  
% Obtained from exceeding the peak threshold  
for k=1:length  
    par3 = floor(index(k)/(m*n))+1;  
    par2 = floor((index(k)-(par3-1)*(m*n))/m)+1;  
    par1 = index(k)-(par3-1)*(m*n)-(par2-1)*m;  
    circleParaXYR = [circleParaXYR;par1,par2,par3];  
    hough_circle(par1,par2)= true; %At this time a lot of centers and radii are gathered. The center of the circle of different circles gathers many points. This is because the given circle is not a standard circle   
end  
  
%The points at the center of each circle are taken as the average, and the exact center and radius of each circle are obtained! 
while size(circleParaXYR,1) >= 1  
    num=1;  
    XYR=[];  
    temp1=circleParaXYR(1,1);  
    temp2=circleParaXYR(1,2);  
    temp3=circleParaXYR(1,3);  
    c1=temp1;  
    c2=temp2;  
    c3=temp3;  
    temp3= r_min+(temp3-1)*step_r;  
   if size(circleParaXYR,1)>1       
     for k=2:size(circleParaXYR,1)  
      if (circleParaXYR(k,1)-temp1)^2+(circleParaXYR(k,2)-temp2)^2 > temp3^2  
         XYR=[XYR;circleParaXYR(k,1),circleParaXYR(k,2),circleParaXYR(k,3)];  %Save the center and radius of the remaining circle 
      else    
      c1=c1+circleParaXYR(k,1);  
      c2=c2+circleParaXYR(k,2);  
      c3=c3+circleParaXYR(k,3);  
      num=num+1;  
      end   
    end  
   end   
      %fprintf(1,'sum %d %d radius %d\n',c1,c2,r_min+(c3-1)*step_r);  
      c1=round(c1/num);  
      c2=round(c2/num);  
      c3=round(c3/num);  
      c3=r_min+(c3-1)*step_r;  
      %fprintf(1,'num=%d\n',num)  
      %fprintf(1,'Center %d %d radius %d\n',c1,c2,c3);     
      para=[para;c1,c2,c3]; %Save the center and radius of each circle  
      circleParaXYR=XYR;  
end  
