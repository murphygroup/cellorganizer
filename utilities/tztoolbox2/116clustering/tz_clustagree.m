function summary = tz_clustagree(cluster1,cluster2,param)
%TZ_CLUSTAGREE The agreement between two 2 clusters.
%   SUMMARY = TZ_CLUSTAGREE(CLUSTER1,CLUSTER2,PARAM) returns the summary of
%   the agreement bewteen two clusters CLUSTER1 and CLUSTER2. PARAM has the
%   following fields:
%       'method' - Evaluation method. It has the following values:
%           'linreg' - linear regression
%               both CLUSTER1 and CLUSTER2 must be [multiple label matrix]
%               or [fuzzy label matrix]. CLUSTER1 will be taken as target
%               variable.
%           'linreg2' - linear regression
%               CLUSTER1 must be a [label vector] and CLUSTER2 must be 
%               [multiple label matrix] or [fuzzy label matrix].             
%           'kappa' - kappa
%               both CLUSTER1 and CLUSTER2 should be a [label vector]
%           'logreg' - logistic regression
%               CLUSTER1 must be a [multiple label matrix] and CLUSTER2 must 
%               be a [multiple label matrix] or [fuzzy label matrix]. CLUSTER1 
%               will be taken as target variable.
%           'mutual' - mutual information
%               both CLUSTER1 and CLUSTER2 should be a [label vector]
%
%   See also

%   17-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 3
    error('Exactly 3 arguments are required')
end


switch param.method
    case 'linreg'
        x = cluster2;
        y = cluster1;
        x(:,sum(x,1)==0) = [];
        [B,res]=tz_transformfeats(x,y,'ln');
        ey = [ones(size(x,1),1) x]*B;
        [tmp,elabel] = max(ey,[],2);
        ey = ml_label2post(elabel,size(ey,2));
        summary.error = sum(res.^2);
        summary.nagree = sum((y(:)==ey(:)).*y(:));
    case 'linreg2'
        x = cluster2;
        y = ml_label2post(cluster1);
        x(:,sum(x,1)==0) = [];
        [B,res]=tz_transformfeats(x,y,'ln');
        ey = [ones(size(x,1),1) x]*B;
        [tmp,elabel] = max(ey,[],2);
        tlabel = cluster1;
        [summary.ncm,summary.pcm,summary.avgacc] = ...
            ml_summaryclassif(tlabel,elabel);
        summary.nagree = sum(diag(summary.ncm));
    case 'kappa'
        [summary.kappa,summary.sd,summary.raw] = xc_kappa(cluster1,cluster2);
    case 'logreg'
        x = cluster2;
        x(:,sum(x,1)==0) = [];
        for i=1:size(cluster1,2)
            y = cluster1(:,i);
            beta = ml_logreg(x,y);
            ps(:,i) = ml_evallogistic(x,beta);
        end
        [tmp,elabel] = max(ps,[],2);
        ey = ml_label2post(elabel,size(ps,2));
        summary.nagree = sum((cluster1(:)==ey(:)).*cluster1(:));
    case 'mutual'
        [summary.info,summary.ncm] = ml_clstmtinfo(cluster1,cluster2);
end
