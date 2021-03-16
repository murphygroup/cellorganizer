function [centers, radii]=DrawBigCell(I,rmin,rmax)
if size(I,3)==3
    I=rgb2gray(I);
end
I=double(I);I=I.^0.4;
I=(I-min(min(I)))/(max(max(I))-min(min(I)));
[centers, radii] = imfindcircles(I,[rmin rmax],'Sensitivity',0.94 );  %The greater the value, the more rounds that are detected
% hold on;viscircles(centers, radii,'EdgeColor','b');