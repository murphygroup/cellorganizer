function [img, proj_orig, eig, locs] = showShapeSpace(model, labels, skipmissing, proj_orig, cm, traces)
    error('This function is depricated. Please use ''showShapeSpaceFigure.m''');

    nimgs = model.cellShapeModel.numimgs;
    cellnums = 1:nimgs;

    if ~exist('labels', 'var') || isempty(labels)
        labels = 1:nimgs;
    end
    
    if ~exist('skipmissing', 'var')
        skipmissing = true;
    end
    
    if ~exist('cm', 'var')
        cm = @jet;
    end
        
    cellimg = model.cellShapeModel.imfunc(2);
    imsize = size(cellimg);

    ndims = 2;    
    gridsize = 3000;

    
     if ~exist('proj_orig', 'var') || isempty(proj_orig)
        if skipmissing
            %remove rows fully populated with nans
            d = model.cellShapeModel.distances_incomplete;
            d_nan = d;
            d_nan(logical(eye(size(d)))) = nan;
            
            rminds = all(isnan(d_nan),1);
            d(rminds,:) = [];
            d(:,rminds) = [];
            
            keepinds = find(~rminds);
            
            labels = labels(keepinds);
            cellnums = cellnums(keepinds);
            
            %remove any other entries with nans
            rminds = any(isnan(triu(d)),1);
            d(rminds,:) = [];
            d(:,rminds) = [];
            
            keepinds = find(~rminds);
            
            labels = labels(keepinds);
            cellnums = cellnums(keepinds);
            
            
            nimgs = size(d,1);
        else
            d = model.cellShapeModel.distances;
        end
        [proj_orig, eig] = cmdscale(d);
     else
         eig = []
     end
        
    
    [ulabels, ~, labelinds ] = unique(labels);
    colors = cm(length(ulabels))*0.8;
    

    
    proj = proj_orig(:,2:-1:1);
   
%     x = proj(:,1);
%     y = proj(:,2);
%     x_ = x ./ (1 + (1 + x.^2 + y.^2).^(1/2));
%     y_ = y ./ (1 + (1 + x.^2 + y.^2).^(1/2));
%     proj(:,1:2) = [x_, y_];
%     
%     [theta, rho] = cart2pol(proj(:,1), proj(:,2));
%     [x,y] = pol2cart(theta, sqrt(rho));
%     proj(:,1:2) = [x,y];
%     
    
%     proj(:,2) = -proj(:,2);
%     
    proj_scale = proj - repmat(min(proj),[nimgs,1]);
    proj_scale = proj_scale/max(proj_scale(:));
    
    locs = round(proj_scale*gridsize);
    
    locs = (locs + max(imsize))./(gridsize+max(imsize)*2) .* gridsize;
    
    if ndims == 2
        locs = [locs ones(nimgs,1)*2];
    end
    
    img = uint8(zeros([repmat(gridsize,[1,ndims]),3]));
    
    
    if exist('traces', 'var')
        traces = traces(all(ismember(traces, keepinds),2),:);
        
        tracemap = ones(max(keepinds), 1) * -1;
        tracemap(keepinds) = 1:length(keepinds);
        
        
        shapeInserter = vision.ShapeInserter('Shape','Lines','BorderColor','Custom', 'CustomBorderColor', uint8([125 125 125]));
        
       lines = int32([locs(tracemap(traces(:,1)),2:-1:1) locs(tracemap(traces(:,2)),2:-1:1)]);
       
       lines = lines + int32(repmat(imsize./4, [size(traces,1) 2]));
%         line = int32([x(i),y(i),x(j), y(j)]);
        img = step(shapeInserter, img, lines);
    end

    for imnum = 1:length(cellnums)
        imnum
        try
            i = cellnums(imnum);
            cellimg = model.cellShapeModel.imfunc(i);
        
            if isempty(cellimg) %|| any(cellimg(:) .* cellbounds(:))
                continue;
            end

            cellmask = cellimg > 0;
            nucmask = cellimg > 1;

            if ndims == 2
                cellmask = sum(cellmask,3) > 0;
                cellmask(cellmask > 0) = 1;
                cellmask = double(cellmask);

                nucmask = sum(nucmask,3) > 0;

                cellmask(nucmask > 0) = 2;

                cellimg = cellmask;
    %             cellimg = sum(cellimg,3);
            end

            cellimg = repmat(cellimg,[1,1,3]);
            cellimg(:,:,1) = cellimg(:,:,1).*colors(labelinds(imnum),1);
            cellimg(:,:,2) = cellimg(:,:,2).*colors(labelinds(imnum),2);
            cellimg(:,:,3) = cellimg(:,:,3).*colors(labelinds(imnum),3);
            
            cellimg = cellimg.*(255/2);
            
            
            idx = find(repmat(sum(cellimg,3), [1,1,size(cellimg,3)])>0);
            [X,Y,Z] = ind2sub(size(cellimg),idx);
%             obj = [X,Y,Z,cellimg(idx)];

% 
%             X = X + locs(i,1);
%             Y = Y + locs(i,1);
%             
%             ind = sub2ind(size(img), X,Y,Z);
%             
%             img(ind) = cellimg(idx);
%             
            X = round(X - mean(X)/2 + locs(imnum,1));
            Y = round(Y - mean(Y)/2 + locs(imnum,2));
            ind = sub2ind(size(img), X, Y, Z);
            img(ind) = cellimg(idx);


% 
%             img = ml_imaddobj2(img,obj,...
%             struct('method','replace','pos',[locs(i,:)]));

        catch
            disp('adsf')
        end
    end
    
    img = flipdim(img,1);
    
    locs = locs(:,2:-1:1);
    locs(:,2) = gridsize - locs(:,2);
end