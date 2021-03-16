function polcodes=tz_cellpolcode(celledge,ra,nucedge,objects,isbin)
%TZ_CELLPOLCODE Coding pixels in polar coordinate system.
%   POLCODES = TZ_CELLPOLCODE(CELLEDGE,RA,NUCEDGE,OBJECTS,ISBIN) returns a
%   matrix of coded pixels. Each row has the form 
%   [normalized distance,angle,gray level] for each pixel. CELLEDGE,
%   NUCEDGE, OBJECTS and ISBIN have the same meaning as those in 
%   TZ_CELLDISTCODE. RA is the angle of the major axis of the cell.
%
%   See also TZ_CELLDISTCODE

%   22-OCT-2004 T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 5
    error('Exactly 5 arguments are required')
end

if any(size(celledge)~=size(nucedge))
    error('sizes of edge images do not match');
end


cellimg=imfill(celledge,'hole');
nucimg=imfill(nucedge,'hole');
objimg=tz_obj2img(objects,size(celledge));

if isbin==1
    objimg=objimg>0;
end

codemask=double(cellimg)-double(nucimg)-double(celledge);

celldist=bwdist(celledge);
nucdist=bwdist(nucedge);
normdistimg=nucdist./(celldist+nucdist);

[I,J]=find(codemask==1);
ind=sub2ind(size(codemask),I,J);

[x,y]=find(nucedge==1);
ca=tz_cellangle([I,J],mean([x,y],1),ra);

polcodes=[nucdist(ind),ca,objimg(ind)];
