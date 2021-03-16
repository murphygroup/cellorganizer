function normdist=tz_cellnormdist(dist1,dist2,dist3,option)
%TZ_CELLNORMDIST Calculate normalized cell distance.
%   NORMDIST = TZ_CELLNORMDIST(DIST1,DIST2,DIST3,OPTION) returns the value
%   of normalized distance upon three distance values:
%       DIST1 - unsigned distance to cell edge
%       DIST2 - signed pixel distance to dna edge
%       DIST3 - signed pixel distance to dna cof  
%   OPTION specifies how to do normalization:
%       'cn' - dnadist/(celldist+dnadist)
%       'mc' - dnadist/(celldist+dnadist) if dnadist>=0
%            - dnadist/(dnadist+dnacofdist) if dandist<0

%   04-Oct-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('Exactly 4 arguments are required')
end

if strcmp(option,'cn')
    normdist=dist2./(dist1+dist2);
    return;
end

if strcmp(option,'mc')
    normdist=dist2;
    ri=find(dist2>=0);
    normdist(ri)=dist2(ri)./(dist1(ri)+dist2(ri));
    li=find(dist2<0);
    normdist(li)=dist2(li)./(abs(dist2(li))+abs(dist3(li)));
else
    error('invalid option');
end