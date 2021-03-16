function OutImg = Normalize0_255(InImg)  
% The pixel value of the superimposed image is normalized by 0-255
ymax=255;ymin=0;  
xmax = max(max(InImg)); %Obtain the maximum value in InImg 
xmin = min(min(InImg)); %Obtain the minimum value in InImg 
OutImg = round((ymax-ymin)*(InImg-xmin)/(xmax-xmin) + ymin); %Normalization and rounding  
end  