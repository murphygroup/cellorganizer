function GfpOutPut=GfpDeal( Path )
FilePath=strcat(Path,'GFP/');
list=dir([FilePath '*.tif']); 
GfpOutPut=cell(size(list,1)-2,1);
for k=3:size(list,1)
 
    FileName=strcat(FilePath,list(k).name);
    Img_Cell=imfinfo(FileName);
    Img_Num=length(Img_Cell);   
    OutImg=cell(Img_Num,1);
    for i=1:Img_Num
        I=imread(FileName,i);  
        originalimg=im2double(I);
        imgtemp = originalimg./max(max(originalimg)); 
        img=255*(imgtemp-min(min(imgtemp)))/(1-min(min(imgtemp)));
        img=uint8(img);
        OutImg{i,1}=img;
    end
    GfpOutPut{k-2,1}=imadd_img(OutImg);
end  


