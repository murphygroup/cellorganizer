function [pt, distance] = findClosest( img, pts )

indices = findn( img~=0 );
point = pts(end,:,:);
point = repmat( point, size(indices,1), 1 );

%d=sqrt((point(1:end,1)-indices(1:end,1)).^2 + ...
%       (point(1:end,2)-indices(1:end,2)).^2 + ...
%       (point(1:end,3)-indices(1:end,3)).^2);

d=sqrt((point(1:end,1)-indices(1:end,1)).^2 + ...
       (point(1:end,2)-indices(1:end,2)).^2 );

pt = indices( find( d==min(d) ), : );
distance = min(d);

if size( pt, 1 ) > 1
  pt = pt(1,:,:);
end
