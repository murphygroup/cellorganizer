function Superposition = Diejia( FilePath1 )
FilePath=strcat(FilePath1,'GFPDeal/');
list=dir(FilePath);  
Superposition=cell(size(list,1)-2,1);
for k=3:size(list,1)  
    filename=strcat(FilePath,list(k).name,'/');
    OutImg=imadd_img(filename);
    Superposition{k-2,1}=OutImg;
end  
end

