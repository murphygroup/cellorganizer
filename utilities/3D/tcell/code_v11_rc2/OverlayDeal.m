function  OverlayOutPut=OverlayDeal( Path )
Img_Cell=imfinfo(strcat(Path,'Overlay.tif'));
Img_Num=length(Img_Cell);
OverlayOutPut=cell(Img_Num,1);
for i=1:Img_Num
    I=imread(strcat(Path,'Overlay.tif'),i);
    OverlayOutPut{i,1}=I;
end

