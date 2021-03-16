function [img,imgfiles]=tz_loadimage(imgdir,ext,loadthre)

%TZ_LOADIMAGE Obsolete. See ML_LOADIMAGE.

error(tz_genmsg('of','tz_loadimage','ml_loadimage'));

if nargin<3
    loadthre=[];
end

imgfiles=ml_dir([imgdir '/*.' ext]);
pos=1;

for i=1:length(imgfiles)
    filenum(i)=tz_getfilenum(imgfiles{i});
end
[sorted,num]=sort(filenum);

if length(sorted)>0
    while any(sorted(1:end-1)-sorted(2:end)==0) & all(num>0)
        pos=pos+1;
        for i=1:length(imgfiles)
            filenum(i)=tz_getfilenum(imgfiles{i},pos);
        end
        if any(filenum<0)
            break
        else
            [sorted,num]=sort(filenum);
        end
    end 
end

imgfiles=imgfiles(num);

for i=1:length(imgfiles)
    img(:,:,i)=ml_readimage([imgdir '/' imgfiles{i}]);
end

if ~isempty(loadthre)
    img(find(img > loadthre)) = 0;
end