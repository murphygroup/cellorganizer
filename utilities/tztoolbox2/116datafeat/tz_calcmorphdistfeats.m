function feat=tz_calcmorphdistfeats(dists,weights,featset)
%TZ_CALCMORPHDISTFEATS Calcuate morphed distance features for a cell.
%   FEAT = TZ_CALCMORPHDISTFEATS(DISTS,WEIGHTS) returns a vector of
%   calculated features of the vector DISTS, which is supposed to be a
%   vector of distances. WEIGHTS is a vector of weights on the distances.
%   DISTAS and WEIGHTS must have the same length. The features are
%   totalweights, mean, variance and skewness.
%   
%   FEAT = TZ_CALCMORPHDISTFEATS(DISTS,WEIGHTS,FEATSET) is under
%   construction.

%   ??-???-2004 Initial write T. Zhao
%   30-OCT-2004 Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('2 or 3 arguments are required')
end

if ~exist('featset','var')
    featset=[];
end

d1=dists(dists<=0);
d2=dists(dists>0);
    
if ~isempty(weights)
    w1=weights(dists<=0);
    w2=weights(dists>0);
else
    w1=[];
    w2=[];
end

feat=[getdistfeats(d1,w1),getdistfeats(d2,w2)]

function feat=getdistfeats(dists,weights)

%function feat=getdistfeats(dists,weights)

if isempty(weights)
    feat(1)=length(dists);
    if isempty(dists)
        feat=[feat,0,0,0];
    else
        v=var(dists);
        feat=[feat,mean(dists),v];
        if v==0
            feat=[feat,0];
        else
            feat=[feat,skewness(dists)];
        end
    end
else
    sw=sum(weights);
    feat(1)=sw;
    mu = sum(dists.*weights)/sw;
    v=sum(((dists - mu).^2).*weights)/sw;
    feat=[feat,mu,v];
    if v==0
        feat=[feat,0];
    else
        t=sum(((dists - mu).^3).*weights)/sw;
        sk=t/v^(3/2);
        feat=[feat,sk];
    end
end
