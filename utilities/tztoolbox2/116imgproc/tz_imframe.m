function img2 = tz_imframe(img,offset,fsize,padvalue)
%TZ_IMFRAME Pad or crop image.
%   IMG2 = TZ_IMFRAME(IMG,OFFSET,FSIZE) pads or crops image at the croners
%   and returns the final image with size FSIZE, which is a length 2 
%   vector. But negative values in fsize will be ignored.
%   OFFSET is the offset of the left top corner. PADVALUE is only useful when
%   the image is expanding in any direction. The default value is 0;

%   13-Sep-2005 Initial write T. Zhao

if nargin < 2
    error('Exactly 2 arguments are required')
end

if ~exist('padvalue','var')
    padvalue = 0;
end

if size(img,3)>1
    for i=1:size(img,3)
        img2(:,:,i) = tz_imframe(squeeze(img(:,:,i)),offset,fsize,padvalue);
    end
    return;
end

if any(fsize) == 0
    img2 =[];
end

if offset(1) > 0
    img = [zeros(offset(1),size(img,2))+padvalue; img];
end

if offset(1) < 0
    img(1:-offset(1),:) = [];
end

if offset(2) > 0
    img = [zeros(size(img,1),offset(2))+padvalue, img];
end

if offset(2) < 0
    img(:,1:-offset(2)) = [];
end
    
sizeDiff = fsize - size(img);

if(fsize(1) > 0)
    if sizeDiff(1) > 0
        img = [img; zeros(sizeDiff(1),size(img,2))+padvalue];
    end
    if sizeDiff(1) < 0
        img(end+sizeDiff(1)+1:end,:) = [];
    end
end

if(fsize(2) > 0)
    if sizeDiff(2) > 0
        img = [img,zeros(size(img,1),sizeDiff(2))+padvalue];
    end

    if sizeDiff(2) < 0
        img(:,end+sizeDiff(2)+1:end) = [];
    end
end

img2 = img;
