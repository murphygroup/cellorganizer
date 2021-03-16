function wfeat = features_textureMR( prot, GLEVELS,dbtype, NLEVELS)

if 1==1
GLEVELS = GLEVELS - 1;

CA = double(prot);

FPL = 27*2;
wfeat = zeros( 1, FPL*NLEVELS);
for k = NLEVELS -1 : -1 : 0;
    [CA,chd,cvd,cdd] = dwt2(CA,'db4');

    A = chd - min(chd(:));
    A = uint8(round(GLEVELS*A/max(A(:))));
    hfeat = ml_texture( A);
    hfeat = [mean(hfeat(1:13,[1 3]),2); mean(hfeat(1:13,[2 4]),2)]';
 
    A = cvd - min(cvd(:));
    A = uint8(round(GLEVELS*A/max(A(:))));
    vfeat = ml_texture( A);
    vfeat = [mean(vfeat(1:13,[1 3]),2); mean(vfeat(1:13,[2 4]),2)]';
 
    hvfeat = (hfeat+vfeat)/2;

    A = cdd - min(cdd(:));
    A = uint8(round(GLEVELS*A/max(A(:))));
    dfeat = ml_texture( A);
    dfeat = [mean(dfeat(1:13,[1 3]),2); mean(dfeat(1:13,[2 4]),2)]';
 
    wfeat(1,(FPL*k)+1:(k+1)*FPL) = [...
        hvfeat dfeat ...
        (sqrt(sum(sum(chd.^2))) + sqrt(sum(sum(cvd.^2))))/2 ...
        sqrt(sum(sum(cdd.^2))) ...
    ];
end

A = CA - min(CA(:));
A = uint8(round(GLEVELS*A/max(A(:))));
dsfeat = ml_texture( A);
dsfeat = [mean(dsfeat(1:13,[1 3]),2); mean(dsfeat(1:13,[2 4]),2)]';

wfeat = [dsfeat wfeat];

return
end

%%%%%%%%%%%%%%%%%%%

if 1==2
tic;

GLEVELS = GLEVELS - 1;

prot = double(prot);

[C,S] = wavedec2(prot,NLEVELS,dbtype);

dsk = appcoef2(C,S,dbtype,NLEVELS);
A = dsk - min(dsk(:));
A = uint8(round(GLEVELS*A/max(A(:))));
dsfeat = ml_texture( A);
dsfeat = [mean(dsfeat(1:13,[1 3]),2); mean(dsfeat(1:13,[2 4]),2)]';

FPL = 27*2;
wfeat = zeros( 1, FPL*NLEVELS);
for k = 0 : NLEVELS-1
    [chd,cvd,cdd] = detcoef2('all',C,S,(NLEVELS-k));

    A = chd - min(chd(:));
    A = uint8(round(GLEVELS*A/max(A(:))));
    hfeat = ml_texture( A);
    hfeat = [mean(hfeat(1:13,[1 3]),2); mean(hfeat(1:13,[2 4]),2)]';
 
    A = cvd - min(cvd(:));
    A = uint8(round(GLEVELS*A/max(A(:))));
    vfeat = ml_texture( A);
    vfeat = [mean(vfeat(1:13,[1 3]),2); mean(vfeat(1:13,[2 4]),2)]';
 
    hvfeat = (hfeat+vfeat)/2;

    A = cdd - min(cdd(:));
    A = uint8(round(GLEVELS*A/max(A(:))));
    dfeat = ml_texture( A);
    dfeat = [mean(dfeat(1:13,[1 3]),2); mean(dfeat(1:13,[2 4]),2)]';
 
    wfeat(1,(FPL*k)+1:(k+1)*FPL) = [...
        hvfeat dfeat ...
        (sqrt(sum(sum(chd.^2))) + sqrt(sum(sum(cvd.^2))))/2 ...
        sqrt(sum(sum(cdd.^2))) ...
    ];

%     wfeat(1,(FPL*k)+1:(k+1)*FPL) = [...
%         hfeat vfeat dfeat ...
%         sqrt(sum(sum(chd.^2))) ...
%         sqrt(sum(sum(cvd.^2))) ...
%         sqrt(sum(sum(cdd.^2))) ...
%     ];
end

wfeat = [dsfeat wfeat];

toc
return
 
end
