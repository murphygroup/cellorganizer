function feat_vec = radInt(image2, imgcent_coordinate)

% The intensity starts only away from the centrosome.
% I will take only the first five.

size_x = size(image2,1);
size_y = size(image2,2);
size_z = size(image2,3);


[cX,cY,cZ] = meshgrid(1:size_x,1:size_y,1:size_z);
eucdist = sqrt(sum(([cX(:)- imgcent_coordinate(2),cY(:)- imgcent_coordinate(1),cZ(:)- imgcent_coordinate(3)]).^2,2));

ranges = ceil(size(image2,3):(max(eucdist)-size(image2,3))/10:(max(eucdist)+size(image2,3)));

for I = 1:size(ranges,2)-1 % The default value is -1 instead of -6. Chose 6 because I want it to be tight.
	rIndices = find((eucdist>ranges(I)).*(eucdist <= ranges(I+1)));
	feat_vec2(I) = sum(image2(rIndices))/size(rIndices,1);
end

feat_vec = feat_vec2(1:5);
