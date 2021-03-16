function [objnum,dists] = tz_genosmparam(osm)
%TZ_GENOSMPARAM Generate parameters for object sampling.
%   OBJNUM = TZ_GENOSMPARAM(OSM) returns a vector of object numbers in 
%   the available clusters. The vector is generated from the object
%   sample model OSM.
%   
%   [OBJNUM,DISTS] = TZ_GENOSMPARAM(...) also returns the generated 
%   distances for the objects.
%   
%   See also TZ_TRAINOSM

%   23-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

objnum=-1;
while any(objnum<0)
    objnum=round(mvnrnd(osm.nmean,osm.nvar));
end

dists = [];

for i=1:length(objnum)
    if objnum(i)==0
        dist = [];
    else
        curdistpara = osm.distpara(osm.occurclst(i),:);
        dist = -1;
        while any(dist<0)
            dist = normrnd(curdistpara(1),curdistpara(2),1,objnum(i));
        end
    end
    dists=[dists dist];
end