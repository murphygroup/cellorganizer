function [prots, nucs] = masks2cropimgs( masks, prot, nuc)
% function [cropimgs,regions] = masks2cropimgs( masks, prot)

bwl = bwlabel(masks==max(masks(:)),4);
s = size(masks);
ma = max(bwl(:));
prots = [];
nucs = [];
regions = [];

for i=1:ma
    [r c] = find(bwl==i);
    minr = min(r);  maxr = max(r);  minc = min(c);  maxc = max(c);
    newIm = uint8(zeros(maxr-minr+1,maxc-minc+1));
    idx = (c-1)*s(1)+r;
    idx2 = (c-min(c))*(max(r)-minr+1)+r-min(r)+1;
    newIm(idx2) = prot(idx);
    
    prots{i} = newIm;

    newIm(idx2) = nuc(idx);
    nucs{i} = newIm;    
end
