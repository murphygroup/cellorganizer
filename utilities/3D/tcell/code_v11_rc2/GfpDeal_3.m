function GfpOutPut=GfpDeal( Path )
FilePath = [Path,'GFP/'];
dir_info = dir(FilePath); 
dir_info(1 : 2) = [];
GfpOutPut=cell(size(dir_info,1),1);
for k=1:size(dir_info, 1)
 
    FileName=strcat(FilePath,dir_info(k).name);
    Img_Cell=imfinfo(FileName);
    Img_Num=length(Img_Cell);   
    OutImg=cell(Img_Num,1);
    for i=1:Img_Num
        I=imread(FileName,i);  
        originalimg=im2double(I);
        imgtemp = originalimg./max(max(originalimg)); 
        img=255*(imgtemp-min(min(imgtemp)))/(1-min(min(imgtemp)));
        img=uint8(img);
        % OutImg{i,1}=img;
        OutImg{i,1}=I;
    end
    constrast_mat = cellfun(@(x) std(double(x(:))), OutImg);
    % center_ind = floor((numel(OutImg) + 1) / 2);
    [~, center_ind] = max(constrast_mat);
    % GfpOutPut{k-2,1}=imadd_img(OutImg(center_ind-3:center_ind+3, 1));
    GfpOutPut{k, 1}=mean(double(cat(3, OutImg{center_ind-5 : center_ind+5})), 3);
    GfpOutPut{k, 1}=mean(double(cat(3, OutImg{:})), 3);
end  


