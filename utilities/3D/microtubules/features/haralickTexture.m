function [ hfeats ] = haralickTexture( imgs, distances )
%cross-channel haralic texture featsures
%
%imgs is a cell array of n x m 2D or 3D images
%
%gj@andrew.cmu.edu
%6/1/14

if isnumeric(imgs)
    imgs = {imgs};
end

%imgs is a cell array of equal sized arrays
nimgs = length(imgs);

h_feats = cell(1, nimgs + ((nimgs^2-nimgs) / 2));

counter = 1;

imgdims = ndims(imgs{1});

for i = 1:nimgs
    for j = i:nimgs
%         disp(num2str([i,j]))
        [~, cr] = cropImg(imgs{i} + imgs{j});
        if imgdims == 3
            h_feats{counter} = cooc3d(imgs{i}(cr(1):cr(2),cr(3):cr(4),cr(5):cr(6)), imgs{j}(cr(1):cr(2),cr(3):cr(4),cr(5):cr(6)), distances);
        elseif imgdims == 2
            h_feats{counter} = cooc3d(imgs{i}(cr(1):cr(2),cr(3):cr(4)), imgs{j}(cr(1):cr(2),cr(3):cr(4)), distances);
        else
            error('images must be 2D or 3D')
        end
        
        counter = counter+1;
    end
end

hfeats = [h_feats{:}];

end

