function [Center,CircleR,Circle_LunKuo,Circle_Center_X,Circle_Center_Y]=BigCellTracking2(DIC,j,Center,CircleR,Circle_Center_X,Circle_Center_Y,width,Superposition,Processdir, param)

if param.verbose
    Processdir_big=strcat(Processdir,'/Big/',num2str(j));
    if ~exist(Processdir_big, 'dir')
        mkdir(Processdir_big);
    end
end
% image1 = imread([file_path1,strcat(int2str(j),'.png')],'png');
image1=DIC;
% image1=histeq(image1);
% Dig out a square area that contains large cells
% img=imcrop(image1,[Circle_Center_X-width/2,Circle_Center_Y-width/2,width,width]);
[img,Position]=Imcrop(image1,Circle_Center_X,Circle_Center_Y,width);
if param.verbose
    imwrite(img,strcat(Processdir_big,'/','01Imcrop.png'),'png');
end
% figure(99),subplot(121),imshow(img)
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
    
    if param.verbose
        figure(99),clf,imshow(img),hold on
        rectangle('Position',[Circle_Center_X-CircleR,Circle_Center_Y-CircleR,2*CircleR,2*CircleR],'Curvature',[1,1]),axis equal
        saveas(gcf,[strcat(Processdir_big,'/'),'02viscircles_circle.png'],'png');
    end
    
    Center=[Circle_Center_X,Circle_Center_Y];
    t=0:0.01*pi:2*pi;
    xxx=cos(t)*CircleR+Circle_Center_X;    
    yyy=sin(t)*CircleR+Circle_Center_Y;
    Circle_LunKuo=[xxx',yyy'];
else
    if param.verbose
        figure(99),clf,imshow(img),hold on
        viscircles(centers, radii,'EdgeColor','b');
        saveas(gcf,[strcat(Processdir_big,'/'),'02viscircles_circle.png'],'png');
    end
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

