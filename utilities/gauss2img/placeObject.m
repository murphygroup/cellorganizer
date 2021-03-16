function [ placedObject ] = placeObject( object, targetSize, targetLocation )
%places object in an matrix of targetSize with the objects centroid at
%targetLocation

%gj 2/15/13

if length(targetSize) == 2
    targetSize = [targetSize 1];
end


if length(size(object)) == 2 & length(targetSize) == 3
    object = padarray(object, [0,0,1]);
end  

%Pad array if it has any singleton dimensions
if sum(size(object) == 1) > 0
    object = padarray(object, double(size(object) == 1));
end

    
objSize = size(object);


if length(targetLocation) == 2
    targetLocation = [targetLocation 1];
end


center = (objSize-1) / 2;

numPx = prod(targetSize);

[y,x,z] = ind2sub(targetSize, 1:numPx);

pxInds = [y', x', z'] + repmat(objSize/2 - targetLocation , [numPx, 1]) ;



if ndims(object) == 2
    a = interp2(object, pxInds(:,1), pxInds(:,2));
else
    a = interp3(object, pxInds(:,2), pxInds(:,1), pxInds(:,3));
end

% a(isnan(a)) = 0;
% 
% 
% placedObject = reshape(a,targetSize);
% 
% 
validinds = ~isnan(a);
placedObject = zeros(targetSize);

placedObject(validinds) = a(validinds);



end

