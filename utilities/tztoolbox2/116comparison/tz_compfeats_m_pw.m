function [pvalues,ts]=tz_compfeats_m_pw(features,method,t)
%TZ_COMPFEATS_M_PW Obsolete.
%
%See also ML_COMPFEATS_M

%ML_COMPFEATS_M_PW pairwise comparison for several matrices
%   [PVALUES,TS]=ML_COMPFEATS_M_PW(FEATS,METHOD,T) performs a multivariate
%   test for all pairs in FEATURES, which is a cell array of matrices.
%   See ml_compfeats_m for details of the parameters METHOD and T.
%   PVALUES is the p-value matirx and TS is the test statistic matrix.

%   T. Zhao

% if ~exist('option','var')
%     option='pw';
% end
% 
% prot_no=length(features);
% 
% for i=1:prot_no
%     for j=(i+1):prot_no
%         [pvalues(j,i),ts(j,i)]=tz_compfeats_m(features{i},features{j},method,t);
%     end
% end

msg = tz_genmsg('of','tz_compfeats_m_pw','ml_compfeats_m_pw');
error(msg);