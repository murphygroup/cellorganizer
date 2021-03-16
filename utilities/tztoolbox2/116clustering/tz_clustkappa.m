function kappa = tz_clustkappa(clustidx1,clustidx2)
%TZ_CLUSTKAPPA Calculate kappa between two clustes.
%   KAPPA = TZ_CLUSTKAPPA(CLUSTIDX1,CLUSIDX2) returns the kappa value of
%   two clusters. CLUSTIDX1 is a vector of cluster labels in the 1st
%   cluster and CLUSTIDX2 is a vector of cluster labels in the 2nd cluster.
%   This function calls XC_KAPPA to calculate kappa value. 
%
%   REFERENCE:
%       X. Chen and R. F. Murphy (2005). Objective Clustering of Proteins 
%       Based on Subcellular Location Patterns. J. Biomed. Biotech. 
%       2005(2):87-95.

%   See also

%   05-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

kappa = xc_kappa(clustidx1,clustidx2);