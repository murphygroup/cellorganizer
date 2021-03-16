function [img, regparam] = cellfile2registeredimg(img, regparam)
    %Center the image so that the center of the nucleus is in the center of
    %the image, the major axis runs left to right, 
    %
    %grj 11/11/13
    %com_align - align the cell such that the center of mass is in the lower left side of the image.
    %
    %grj 8/26/14
    %bottom_align - align the cell so the bottoms align, rather than the COM in the Z-direction
    
    if ~exist('regparam', 'var')
        regparam = struct();
    end
    
    regparam = ml_initparam(regparam, ...
        struct('zflip', true, ...
               'z_align', 'com', ... %can also be 'com', or 'largest_slice'
               'shape_to_model', 'both', ...  %can also be 'cell', or 'nuc' 
               'com_align', 'cell')); %can also be 'cell', or 'nuc' 
           
    if strcmpi(regparam.shape_to_model, 'nuc') && ~strcmpi(regparam.com_align, 'body')
        warning('regparam.shape_to_model set to ''nuc''. Setting regparam.com_align to ''body''')
        regparam.com_align = 'body';
    end
        
    if ischar(img)
        fname = img;
        
        if ~exist(fname, 'file')
            img = [];
            return
        end

        load(fname, 'segdna', 'segcell')

        img = zeros(size(segdna));

        boundimg = segdna;

        %cell boundary is a bunch of 1's
        if exist('segcell', 'var')
            img(segcell>0) = 1; 
            boundimg = segdna + segcell;
        end

        imbounds = true(size(segdna));
        imbounds(2:end-1,2:end-1,:) = 0;

        %nucleus is a bunch of 2's
        if isempty(segdna) | all(~segdna(:)) | any(boundimg(imbounds(:)))
            img = [];
            return
        end

        img(segdna>0) = 2;
    end

    switch regparam.shape_to_model
        case 'cell'
            img = double(img>=1);
        case 'nuc'
            img = double(img>=2);
        otherwise
            %do nothing
    end
    
    img = ml_findmainobj(img);
    
    px_per_slice = squeeze(sum(sum(img>0,1),2));
    %if the image is 3-dims make sure it's not upside-down
    
    if ~isfield(regparam, 'flipdim')

        if strcmpi(regparam.com_align, 'nuc')
            c1 = regionprops(img>0, 'Centroid');
            c2 = regionprops(img>1, 'Centroid');
            regparam.flipdim = (c1.Centroid(1) - c2.Centroid(1)) < 0;
        else
            inds = find(img>0);
            [x,y,~] = ind2sub(size(img), inds);
            skew = skewness([x,y]);

            skew = skew(1:2);

            regparam.flipdim = skew < 0;
        end
    end

    for i = 1:length(regparam.flipdim)
        if regparam.flipdim(i)
            img = flip(img, i);
        end
    end
    
    
    if regparam.zflip 
        if ndims(img) == 3
            stats = regionprops(img>0);
            centroid = stats.Centroid;

            %make sure the image is not upside-down
            slice_inds = find(px_per_slice);
            if (slice_inds(end) + slice_inds(1)) / 2 < centroid(3)
                img = flipdim(img, 3);
            end
        end

    end

    %if an angle of rotation is not defined, find  it
    if ~isfield(regparam, 'angle')
        
%         if strcmpi(regparam.com_align, 'nuc')
%             c1 = regionprops(sum(img,3)>0, 'Centroid');
%             c2 = regionprops(sum(img,3)>1, 'Centroid');
%             
%             %angle the iamge diagonally to save space
%             regparam.angle = 45 + atan2d(c1.Centroid(1) - c2.Centroid(1), c1.Centroid(2) - c2.Centroid(2));
%         else
            angle = regionprops(sum(img,3)>0, 'orientation');
            regparam.angle = 45 + angle.Orientation;
%         end
    end
    
    img = imrotate(img, -regparam.angle);
    img = double(ml_findmainobj(img>0)) + double(ml_findmainobj(img>1));

    %             figure, imshow(sum(imrot,3),[])
    %rotate and crop the image to figure out how big it is
    if ~isfield(regparam, 'croprange')
        [img, regparam.croprange] = cropImg(img, 5);
    else
        img = img(regparam.croprange(1):regparam.croprange(2), regparam.croprange(3):regparam.croprange(4),:);
    end

    imsize = size(img);
    
    if length(imsize) == 2
        imsize = [imsize 1];
    end

    px_per_slice = squeeze(sum(sum(img>0,1),2));
    [~, maxind] = max(px_per_slice);

    %sometimes there are disconnected regions. In practice they are quite
    %small and we ignore them
    if ~isfield(regparam, 'offset')
        switch regparam.com_align
            case 'cell'
                regions = regionprops(img>0);
            case 'nuc'
                regions = regionprops(img>1);
            otherwise
                regions = regionprops(img>0);
        end
                
        [~, objind] = max([regions(:).Area]);

        centroid = round(regions(objind).Centroid );
        
        if size(img,3) == 1;
            centroid(3) = 0.5;
        end
        
        
        switch regparam.z_align
            case 'largest_slice'
                centroid(3) = maxind;
            case 'bottom'
                centroid(3) = find(px_per_slice > 0, 1, 'first');
            case 'com'
                %do nothing
            otherwise
                %do nothing
        end
        centroid(1:2) = centroid(2:-1:1);

        offset = round(imsize - centroid*2);

        if imsize(3) == 1
            offset(3) = 0;
        end

        regparam.offset = offset;
    end
    
    offset = regparam.offset;
    
    for j = 1:length(offset)
        pad = zeros(1, length(imsize));
        pad(j) = abs(offset(j));
        if offset(j) == 0
            continue
        elseif offset(j) > 0
            str = 'pre';
        else
            str = 'post';
        end
        img = padarray(img, pad, str);
    end
    

    
    if isfield(regparam, 'downsample')
        %D. Sullivan 5/5/14 - don't downsample if you don't have a 3rd
        %dimension.
        if size(img,3)==1
            regparam.downsample(3) = 1;
        end
        img = ml_downsize(img, regparam.downsample, 'nearest');
    end
end