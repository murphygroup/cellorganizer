function primitives = getBox3(img,name,resolution,primitives,param)


if iscell(img)
    imgtmp = double(img{1});
    for i = 1:length(img)
        imgtmp = double(img{i})+imgtmp;
    end
    img = imgtmp>0;
elseif ndims(img)~=3
    warning('Unable to make bounding box - Unsupported dimension');
    return
end

if nargin<4 || isempty(primitives)
    primitives = struct;
    primitives.(name).list = [];
end

if param.framework_cropping
    %trim image
    img2d = squeeze(sum(img,3));
    xrange = find(sum(img2d,1)>0);
    yrange = find(sum(img2d,2)>0);
    zrange = find(sum(sum(img,1),2));

    imgcrop = img(min(yrange):max(yrange),min(xrange):max(xrange),min(zrange):max(zrange));
else
    xrange = [1, size(img,1)];
    yrange = [1, size(img,2)];
    zrange = [1, size(img,3)];
end

%calculate size required
if param.framework_cropping
    imgsize =size(imgcrop);
else
    imgsize =size(img);
end
imgsize = [imgsize(2),imgsize(1),imgsize(3)];
% padding = [1,1,1];
padding = [2,2,1];
% padding = [3,3,1];
Boxsize = (imgsize.*resolution)/2+padding;
    

%calculate position

ind = length(primitives.(name).list)+1;
primitives.(name).list(ind).name = name;
primitives.(name).list(ind).type = 'cube';
%D. Sullivan 10/25/14 - no longer need the shift since the padding is now removed from the cell meshing
primitives.(name).list(ind).position = [0,0,0]+[ceil(mean(xrange)),ceil(mean(yrange)),ceil(mean(zrange))].*resolution;%+padding;
primitives.(name).list(ind).rotation = [0,0,0];
primitives.(name).list(ind).scale = Boxsize;
primitives.(name).list(ind).resolution = [1,1,1];
primitives.(name).ordinal = 0;

%D. Sullivan - make sure the bounding box is also scaled to
%match the scaling done on the CP (ensuring that the CP doesn't
%overlap with the nuc)
param = ml_initparam(param,struct('adjustscale',1.025));
param.adjustsize = 1;
if isfield(param,'adjustsize') && param.adjustsize
    %         tmpverts = FV.vertices.*1.1;
%     tmpscale = primitives.(name).list(ind).scale.*param.adjustscale;
    Boxsize = primitives.(name).list(ind).scale.*param.adjustscale;
    %need to shift the size adjusted vertices to the correct postions.
    %                 FV.vertices = tmpverts-repmat(max(tmpverts-FV.vertices)/2,size(tmpverts,1),1);
    primitives.(name).list(ind).position = primitives.(name).list(ind).position-(Boxsize-primitives.(name).list(ind).scale)/2;
    primitives.(name).list(ind).scale = Boxsize;%Make sure this is the same as in addMeshObjects/2 (remember scale is x2)
    
    %         FV.vertices = FV.vertices+1;
end


%Adding volume and surface area fields
Boxsize = Boxsize.*2;
primitives.(name).totvol =  prod(Boxsize);
primitives.(name).totsa = 2*(prod(Boxsize(1:2))+prod(Boxsize(2:3))+prod([Boxsize(1),Boxsize(3)]));



end

