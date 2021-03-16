function combobjpsize=tz_calcobjpsize(combobjects)
%TZ_CALCOBJPSIZE Calculate number of pixels in each object.
%   COMBOBJPSIZE = TZ_CALCOBJPSIZE(COMBOBJECTS) returns a vector of the
%   sizes of objects in the one-level cell array COMBOBJECTS. 

%   ??-???-???? Initial write T. Zhao
%   03-NOV-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

nobj=length(combobjects);

for i=1:nobj
    combobjpsize(i)=size(combobjects{i},1);
end