function pts2 = tz_expandpts(pts,order)
%TZ_EXPANDPTS Obsolete. See ML_EXPANDPTS.
%   PTS = TZ_EXPANDPTS(PTS,ORDER) returns a matrix, each row of which is
%   extended from a point in PTS by an order of the absolute value of 
%   ORDER. If ORDER is less than 0, zero order will be added.
%   
%   See also

%   07-Apr-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_expandpts','ml_expandpts'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

pts2=[];
for i=1:abs(order)
    pts2=[pts2,pts.^i];
end

if order<=0
    pts2=[pts2,ones(size(pts2,1),1)];
end
