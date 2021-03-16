function objfeats=tz_selobjfeatpmcf2comb(comobjfeats,sel)
%TZ_SELOBJFEATPMCF2COMB Combine selected object features.
%   OBJFEATS = TZ_SELOBJFEATPMCF2COMB(COMOBJFEATS,SEL) returns a feature
%   matrix of objects in the one-level cell array comobjfeats. Only cells
%   with indices in the vector SEL will be selected. The last column in
%   each feature matrix of comobjfeats must be the vector of cell indices.
%   
%   See also

%   24-MAR-2004 Initial write T. Zhao
%   04-NOV-2004 Modified T. Zhao
%       - add comments
%       - change funciton name tz_getobjfeats --> tz_selobjfeatpmcf2comb
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

nclass=length(comobjfeats);
objfeats=[];

for i=1:nclass
    if isempty(sel)
        objfeats=[objfeats;comobjfeats{i}];
    else
        for j=sel
            objfeats=[objfeats;comobjfeats{i}(comobjfeats{i}(:,end)==j,:)];
        end
    end
end

