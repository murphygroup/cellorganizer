function [ImgNew,Position]=Imcrop(Img,Circle_Center_X,Circle_Center_Y,Width)
Position=zeros(1,4);
Position(1)=Circle_Center_X-Width/2;
Position(2)=Circle_Center_Y-Width/2;
Position(3)=Width;
Position(4)=Width;
% Border processing
if Position(1)<1
    Position(1)=1;
end
if Position(1)+Width>size(Img,2)
    Position(1)=size(Img,2)-Width;
end
if Position(2)<1
    Position(2)=1;
end
if Position(2)+Width>size(Img,1)
    Position(2)=size(Img,1)-Width;
end

ImgNew=imcrop(Img,Position);
