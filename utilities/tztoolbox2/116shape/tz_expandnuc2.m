function cellpts = tz_expandnuc2(cellcode,angles,rdist)
%TZ_EXPANDNUC2 Expand a nucleus to generate a cell.
%   CELLPTS = TZ_EXPANDNUC2(CELLCODE,ANGLES,RDIST) returns an array of
%   points that represent the cell expanded from the nucleus described 
%   in CELLCODE ('nuchitpts'). ANGLES is a vector and it specifies the
%   set of angles for transformation. RDIST is also a vector and it must
%   have the same size as ANGLES. Each element in RDIST is the ratio of the
%   distance from nucdist to celldist, where nucdist is the distance from
%   nuclear center to nuclear boundary and celldist is the distance from
%   nuclear center to cell boundary. See TZ_PARSECELL for the sturcture of
%   CELLCODE.
%   
%   See also TZ_EXPANDNUC

%   28-Mar-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

if isempty(angles)
    angles=1:360;
end

dists=cellcode.nucdist(angles);
celldist=dists./rdist;

for i=1:length(dists)
    pts=tz_getlinept2(cellcode.nuccenter(1:2),angles(i),celldist(i),1);
    cellpts(i,:)=pts(end,:);
end
