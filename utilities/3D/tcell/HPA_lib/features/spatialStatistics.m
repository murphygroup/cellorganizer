function stats = spatialStatistics(img1, img2, DSF)

s = size(img1);

mask = zeros(size(img1));
mask(1:DSF:s(1),1:DSF:s(2)) = 1;
[r1,c1] = find(img1 & mask);
[r2,c2] = find(img2 & mask);

coord1 = [r1,c1];
coord2 = [r2,c2];

odist = pdist( [coord1; coord2]);
odist = squareform(odist);

% We do not want to mix distances from same channel with distances from other channel
odist = odist(1:size(coord1,1),size(coord1,1)+1:end);

odist = triu(odist+1, 1) - 1;
odist(odist==-1) = [];

NBINS = 39;
bsize = 1/NBINS;

odist = odist - min(odist);
odist = odist / max(odist);

[c b] = hist(odist,0:bsize:1);

stats = c / sum(c) / bsize;

return
