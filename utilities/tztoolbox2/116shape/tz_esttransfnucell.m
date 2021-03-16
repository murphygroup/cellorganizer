function [A,ress] = tz_esttransfnucell(cellcode,option,order)
%TZ_ESTTRANSFNUCELL Estimate transformation from nucleus to cell.
%   A = TZ_ESTTRANSFNUCELL(CELLCODE,OPTION,ORDER) returns the
%   transformation matrix from nucleus to cell. All information of the
%   nucleus and cell is contained in the structure CELLCODE. OPTION
%   specifies transformation method:
%       'af' - affine transformation
%       'ac' - rotation, scaling and translation
%   ORDER is only useful for affine transformation.
%
%   [A,RESS] = TZ_ESTTRANSFNUCELL(...) also returns the transformation
%   error.
%   
%   See also TZ_DEFORMNUCELL

%   23-Apr-2005 Initial write Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University
 

%move center to the original point
pts1=ml_addrow(cellcode.nuchitpts,-cellcode.nuccenter);
pts2=ml_addrow(cellcode.cellhitpts,-cellcode.nuccenter);
cellpts=ml_addrow(cellcode.cellcontour,-cellcode.nuccenter);

%get transformation matrix
switch(option)
case 'af'   %affine transformation
    A=tz_estransf_2d(pts1,pts2,order);
case 'ac' %active model
    [a,t]=tz_alignshape(pts1,pts2);
    A=[a(1) a(2) 0;-a(2) a(1) 0;[t 1]];
end
    
%nucleus transformation
tpts=tz_deformpts_2d(pts1,A);

%move to positive part
mintpts=min(tpts,[],1);
tpts=ml_addrow(tpts,-mintpts)+1;

%generate closed contour
conpts=tz_showpts_2d(round(tpts),'ln',1);

%solid object image of the transformed nucleus
imgsize=max(conpts,[],1);
edgeimg=tz_obj2img(conpts,imgsize);
objimg=imfill(edgeimg,'hole');

%get major angle of the transformed nucleus
tmangle=tz_bwmajorangle(objimg);


%get center of the transformed nucleus
[x,y]=find(objimg==1);
tcenter=mean([x,y],1);
ztnuc=ml_addrow(tpts,-tcenter);
tcenter=tcenter+mintpts-1;

%calculate angles between edge points and major axis
tas=tz_ptlnangle(ztnuc,[0 0 tmangle]);

%closest to major axis
[mintas,mintasidx]=min(abs(tas));

ztcell=ml_addrow(cellpts,-tcenter);
cas=tz_ptlnangle(ztcell,[0 0 tmangle]);

for i=1:length(tas)
    [tmp,caidx(i)]=min(abs(tas(i)-cas));
end

% tpts=tz_normapts(tpts,[]);
[dists,d1,d2] = tz_diffhitpts(ztnuc,ztcell(caidx,:));

if mintasidx>1
    ress=[dists(mintasidx:end);dists(1:mintasidx-1)];
else
    ress=dists;
end

% rdist=dists./(d2-normnucdist');
% dists=rdist.*(d2-normnucdist');