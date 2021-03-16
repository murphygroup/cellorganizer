function corrmats = tz_imcorrtex2(img,mask,bin,option)
%TZ_IMCORRTEX2 Co-occurence matrix of an image.
%   TEX = TZ_IMCORRTEX2(IMG,MASK,BIN) returns a co-occurence matrix 
%   of the image IMG. which will be masked if MASK is not empty. BIN
%   specifies how many bins will be used for the matrix if BIN is not empty
%   OPTION is a string. BIN can also be negative. If it is negative, the
%   image will be binned to -BIN bins.
%   
%   TEX = TZ_IMCORRTEX2(IMG,MASK,BIN,OPTION) specifies additional options.
%   If OPTION is 'rm0', dark pixels (gray level 0) will not be counted.
%   
%   See also

%   22-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('3 or 4 arguments are required')
end

if nargin < 4
    option = '';
end

if ~strcmp(class(img),'double')
    img = double(img);
end
    
if ~isempty(mask)
    maskIndices = find(mask==0);
end

if ~isempty(mask)
    img(maskIndices) = -1;
end

if strcmp(option,'rm0')
    img(img==0) = -1;
end

validIndices = find(img>=0);
if isempty(validIndices)
    if ~isempty(bin)
        tex = zeros(abs(bin));
    else
        tex = [];
    end
    corrmats = tex;
    return
end

minlevel=min(img(validIndices));
maxlevel=max(img(validIndices));

if ~isempty(bin)
    if bin==0
        error('BIN should not be 0');
    end
    if bin<0
        bin = -bin;
        
        if maxlevel == minlevel
            if maxlevel>=bin
                img(validIndices) = bin;
            end
        else
            
            img(validIndices) = round((img(validIndices)-minlevel)/ ...
                (maxlevel-minlevel)*(bin-1));
        end
    else
        img(img>=bin) = -1;
    end
else
    bin = maxlevel+1;
end

img = img+1;

shiftImage(:,:,1) = tz_imtranslate(img,[0,1]);
shiftImage(:,:,2) = tz_imtranslate(img,[1,0]);
shiftImage(:,:,3) = tz_imtranslate(img,[1 1]);
shiftImage(:,:,4) = tz_imtranslate(img,[-1 1]);


pixelPair = [];
for k=1:size(shiftImage,3)
    currentShiftImage = squeeze(shiftImage(:,:,k));
    validIndices = find(img>0 & currentShiftImage>0);
    pixelPair = [img(validIndices) currentShiftImage(validIndices)];

    tex = zeros(bin,bin);
    if ~isempty(pixelPair)
        corrIndices = sub2ind([bin bin],pixelPair(:,1),pixelPair(:,2));
        corrHist = hist(corrIndices,max(corrIndices)-min(corrIndices)+1);    
        tex(unique(corrIndices)) = corrHist(corrHist>0);
        tex = tex+tex'-diag(diag(tex));
    end
    corrmats(:,:,k) = tex;
end