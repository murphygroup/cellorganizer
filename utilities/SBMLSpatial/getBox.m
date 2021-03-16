function primitives = getBox(img,name,resolution,primitives,param)


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

xrange = [1, size(img,2)];
yrange = [1, size(img,1)];
zrange = [1, size(img,3)];

%calculate size required
imgsize =size(img);
imgsize_xyz = [imgsize(2),imgsize(1),imgsize(3)];
% use_padding = false;
use_padding = true;
if use_padding
    % padding = [1,1,1];
    padding = [2,2,1];
    % padding = [3,3,1];
else
    padding = [0,0,0];
end
Boxsize = (imgsize_xyz.*resolution)/2+padding;


%calculate position

ind = length(primitives.(name).list)+1;
primitives.(name).list(ind).name = name;
primitives.(name).list(ind).type = 'cube';
%D. Sullivan 10/25/14 - no longer need the shift since the padding is now removed from the cell meshing
primitives.(name).list(ind).position = [0,0,0]+[ceil(mean(xrange)),ceil(mean(yrange)),ceil(mean(zrange))].*resolution;%+padding;
primitives.(name).list(ind).rotation = [0,0,0];
primitives.(name).list(ind).scale = Boxsize;
primitives.(name).ordinal = 0;


%Adding volume and surface area fields
Boxsize = Boxsize.*2;
primitives.(name).totvol =  prod(Boxsize);
primitives.(name).totsa = 2*(prod(Boxsize(1:2))+prod(Boxsize(2:3))+prod([Boxsize(1),Boxsize(3)]));



end

