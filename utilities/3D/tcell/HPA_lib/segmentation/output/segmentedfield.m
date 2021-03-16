function output = segmentedfield(nuc,cell,regions)

[c b] = imhist(nuc);
[a i] = max(c);
nuc = nuc - b(i);
nuc = (255/max(nuc(:)))*nuc;

[c b] = imhist(cell);
[a i] = max(c);
cell = cell - b(i);
cell = (255/max(cell(:)))*cell;

output(:,:,3) = nuc;
output(:,:,2) = cell;
output(:,:,1) = 32*uint8(regions) + 128*uint8( del2(double(regions))~=0);
