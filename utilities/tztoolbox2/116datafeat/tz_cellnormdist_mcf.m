function normdists=tz_cellnormdist_mcf(celldists,dnadists,option)
%TZ_CELLNORMDIST_MCF Calculate normalized distances for mcf
%   NORMDISTS = TZ_CELLNORMDIST_MCF(CELLDISTS,DNADISTS,OPTION) returns a
%   cell array of normalized distances. CELLDISTS{I}{J}{K} is the distance
%   between the Kth object and the cell edge in Jth cell of the Ith class.
%   DNADIST{I}{J}{K} has two columns. The first column is the distance to 
%   dna edge and the second column is the distance to dna center.
%   See TZ_CELLNORMDIST for more details.

%   11-Sep-2005 Initial write T. Zhao
%   ??-???-???? Initial T. Zhao
%   31-OCT-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

nclass=length(celldists);

for i=1:nclass
    ncell=length(celldists{i});
    
    for j=1:ncell
        [i j]
        nobj=length(celldists{i}{j});
        for k=1:nobj
            normdists{i}{j}{k}=tz_cellnormdist(celldists{i}{j}{k},dnadists{i}{j}{k}(:,1),dnadists{i}{j}{k}(:,2),option);
        end
    end
end