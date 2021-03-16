function img = tp_findmainregion(img)
% Find the largest connected region in image. IMG is a binary matrix.

if ndims(img) > 2
    [L,n] = bwlabeln(img,18);
else
    [L,n] = bwlabel(img,8);
end
for i = 1:n
    Area(i) = nnz(L==i);
end
[y,maxl] = max(Area);
img(L~=maxl) = 0;