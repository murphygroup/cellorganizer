function tz_showobjcom(cellobjcom,sel)
%TZ_SHOWOBJCOM Show histogram like object composition.
%   TZ_SHOWOBJCOM(CELLOBJCOM) show the object composition for each class
%   in CELLOBJCOM, which is a one-level cell array of object compostions
%   in each cell.
%   
%   TZ_SHOWOBJCOM(CELLOBJCOM,SEL) selects cells by the indexing vector 
%   SEL for showing.

%   ??-???-2004 Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required')
end

nclass=length(cellobjcom);
for i=1:nclass
    subplot(nclass,1,i);
    if ~exist('sel','var')
        bar(sum(cellobjcom{i}(:,1:end-1)));
    else
        bar(sum(cellobjcom{i}(sel,1:end-1)));
    end
end