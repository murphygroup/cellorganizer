function [feats, names] = hpatilefeatures( prot, dna, coords)

GLEVELS = 32;
dbtype = 'db4';
NLEVELS = 10;

bw = dna>255*graythresh(dna);

feats = zeros(size(coords,1),1+26+3+(26+2*27*NLEVELS)+18);
for i=1:size(coords,1)
    rs = coords(i,1);
    re = coords(i,2);
    cs = coords(i,3);
    ce = coords(i,4);
    protim = prot(rs:re,cs:ce);
    dnaim  =  dna(rs:re,cs:ce);
    bwim   =   bw(rs:re,cs:ce);

    ffeats = sum(bwim(:));
    ofeats = features_nuclearOverlap( protim, dnaim);
    tfeats = features_texture( protim, GLEVELS);

    wfeats = features_textureMR(protim,GLEVELS,dbtype,NLEVELS);

    [names,tasfeats,slfnames] = ml_tas( protim);
    feats(i,:) = [ffeats ofeats tfeats wfeats tasfeats];
end

names(1) = {'filter'};
names(2:4) = {'overlap'};
names(5:30) = {'texture'};
names(31:size(feats,2)) = {'wavelet'};
