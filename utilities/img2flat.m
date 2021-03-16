function cellimg = img2flat(cellimg, color)
    ndims = 2;

    if isempty(cellimg) %|| any(cellimg(:) .* cellbounds(:))
        return
    end
    
    if ~exist('color', 'var') | isempty(color)
        color = [1, 0, 0; ...
                  0, 0, 1];
    end

    cellmask = cellimg > 0;
    nucmask = cellimg > 1;
    
    if size(color,2) == 2
        cellmask = sum(double(cellimg > 0),3)/size(cellimg,3);
        nucmask = sum(double(cellimg > 1),3)/size(cellimg,3);

        cellmask = repmat(cellmask,[1,1,3]);
        cellmask(:,:,1) = cellmask(:,:,1).*color(1,1);
        cellmask(:,:,2) = cellmask(:,:,2).*color(1,2);
        cellmask(:,:,3) = cellmask(:,:,3).*color(1,3);


        nucmask = repmat(nucmask,[1,1,3]);
        nucmask(:,:,1) = nucmask(:,:,1).*color(2,1);
        nucmask(:,:,2) = nucmask(:,:,2).*color(2,2);
        nucmask(:,:,3) = nucmask(:,:,3).*color(2,3);

        cellmask(repmat(sum(nucmask,3), [1,1,3]) > 0) = 0;

        cellimg = cellmask + nucmask;
    else
    
        cellmask = sum(cellmask,3) > 0;
        cellmask(cellmask > 0) = 1;
        cellmask = double(cellmask);

        nucmask = sum(nucmask,3) > 0;

        cellmask(nucmask > 0) = 2;

        cellimg = cellmask;
%         cellimg = cellmask/2;
%             cellimg = sum(cellimg,3);
    

        cellimg = repmat(cellimg,[1,1,3]);
        cellimg(:,:,1) = cellimg(:,:,1).*color(1);
        cellimg(:,:,2) = cellimg(:,:,2).*color(2);
        cellimg(:,:,3) = cellimg(:,:,3).*color(3);
    end

end