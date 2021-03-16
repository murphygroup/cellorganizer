function [Center,CircleR,Circle_LunKuo,Circle_Center_X,Circle_Center_Y]=BigCellTracking(DIC,Center,CircleR,Circle_Center_X,Circle_Center_Y,width)
image1 = DIC;
[img,Position]=Imcrop(image1,Circle_Center_X,Circle_Center_Y,width);
if size(img,3)==3
    I=rgb2gray(img);
else
    I=img;
end
I=double(I);I=I.^0.5;
I=(I-min(min(I)))/(max(max(I))-min(min(I)));
% subplot(122),imshow(I)
[centers, radii] = imfindcircles(I,[15 35],'Sensitivity',0.95);
centers=round(centers);
radii=round(radii);

% below discuss three cases where a tester may fail:
%1: if the number of circles detected is empty, use the previous round of data
%2: if the detected circle radius is 1.5 times longer than the previous circle radius, the circle data of the previous frame is used
%3: if the distance between the center of the circle and the center of the previous frame is more than 1 times the radius of the circle, the circle data of the previous frame will be used
% Determines whether the circle is detected
if isempty(centers)
    % Not detected, case 1
    max_R=CircleR;
    CircleR=CircleR;
    Circle_Center_X=Circle_Center_X;
    Circle_Center_Y=Circle_Center_Y;
%     figure(99),clf,imshow(img),hold on
%     rectangle('Position',[Circle_Center_X-CircleR,Circle_Center_Y-CircleR,2*CircleR,2*CircleR],'Curvature',[1,1]),axis equal
    Center=[Circle_Center_X,Circle_Center_Y];
    t=0:0.01*pi:2*pi;
    xxx=cos(t)*CircleR+Circle_Center_X;    
    yyy=sin(t)*CircleR+Circle_Center_Y;
    Circle_LunKuo=[xxx',yyy'];
else
%     figure(99),clf,imshow(img),hold on
%     viscircles(centers, radii,'EdgeColor','b');
    circleParaXYR=[centers(:,2),centers(:,1),radii];
    % The center of the next frame is determined by detecting the center of the circle and the center of the previous frame (the center of the center of the previous frame is actually the center of the cut)
    [Circle_Center_X2,Circle_Center_Y2,max_R,min_Distance]=GetCircle_ByCenter_Distance(circleParaXYR,width/2,width/2);
        Circle_Center_X=Circle_Center_X2+Position(1);
        Circle_Center_Y=Circle_Center_Y2+Position(2);
    if (max_R>=1.5*CircleR) |  (min_Distance>=1.5*CircleR)  % If the circle radius of the frame is greater than 1.5 times of the previous frame, the circle data of the previous frame is used
        max_R=CircleR;
        Circle_Center_X=Center(1);
        Circle_Center_Y=Center(2);
    end
    t=0:0.01*pi:2*pi;
    xxx=cos(t)*max_R+Circle_Center_X;    
    yyy=sin(t)*max_R+Circle_Center_Y;
    Center=[Circle_Center_X,Circle_Center_Y];
    CircleR=max_R;
    Circle_LunKuo=[xxx',yyy'];
end

