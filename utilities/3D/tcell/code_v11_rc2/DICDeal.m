function DICOutPut=DICDeal( Path  )
dir_info = dir(Path);
dir_info(1:2) = [];
filenames = {dir_info.name};
file_pattern = '.*DIC.*\.tiff?$';
filename = filenames{~cellfun(@isempty, regexpi(filenames, file_pattern))};
Img_Cell=imfinfo(strcat(Path,filename));
Img_Num=length(Img_Cell);

DICOutPut=cell(Img_Num,1);
for i=1:Img_Num
    I=imread(strcat(Path,filename),i);
    originalimg=im2double(I);
    imgtemp = originalimg./max(max(originalimg)); 
    img=255*(imgtemp-min(min(imgtemp)))/(1-min(min(imgtemp)));
    img=uint8(img);
    DICOutPut{i,1}=img;
end
end
