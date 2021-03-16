function y = tz_evalsvmreg(sample,regmodel)
%TZ_EVALSVMREG Obsolete

%function y = tz_evalsvmreg(sample,regmodel)
%OVERVIEW
%   evaluate svm regerssion
%PARAMETERS
%   sample - testing data
%   regmodel - trained bpnn
%RETURN
%   y - predicted targets
%DESCRIPTION
%   
%HISTORY
%   15-June-2005 Initial write SHANNCC
%SEE ALSO
%   

error(tz_genmsg('of','tz_evalsvmreg','ml_evalsvmreg'));

if ~strcmp(regmodel.modeltype,'svm')
    error('wrong model type: not a support vector machine');
end

if isfield(regmodel,'prep')
    if ~isempty(regmodel.prep.featidx)
        sample=sample(:,regmodel.prep.featidx);
    end
    if regmodel.prep.zscore
        sample=tz_zscore(sample,regmodel.prep.zmean,regmodel.prep.zsdev);
    end
end

% Summarize network performance     Sam testing part, change to fwd(net...)
[y_score, y] = my_fwd(regmodel.trained, sample);
% [nmax, y] = max(y,[],2) ;

if isfield(regmodel,'postp');
    if isfield(regmodel.postp,'ctg')
        if regmodel.postp.ctg==1 %categorize
            [nmax, y] = max(y,[],2) ;
        end
    end
end