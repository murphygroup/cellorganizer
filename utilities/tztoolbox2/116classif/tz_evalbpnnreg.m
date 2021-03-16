function y = tz_evalbpnnreg(sample,regmodel)
%TZ_EVALBPNNREG Obsolete

%function y = tz_evalbpnnreg(sample,regmodel)
%OVERVIEW
%   evaluate bpnn regerssion
%PARAMETERS
%   sample - testing data
%   regmodel - trained bpnn
%RETURN
%   y - predicted targets
%DESCRIPTION
%   
%HISTORY
%   14-May-2005 Initial write TINGZ
%SEE ALSO
%   
error(tz_genmsg('of','tz_evalbpnnreg','ml_evalbpnnreg'))

if ~strcmp(regmodel.modeltype,'network')
    error('wrong model type: not a network');
end

if isfield(regmodel,'prep')
    if ~isempty(regmodel.prep.featidx)
        sample=sample(:,regmodel.prep.featidx);
    end
    if regmodel.prep.zscore
        sample=tz_zscore(sample,regmodel.prep.zmean,regmodel.prep.zsdev);
    end
end

% Summarize network performance
y = mlpfwd(regmodel.trained, sample);

if isfield(regmodel,'postp');
    if regmodel.postp.ctg %categorize
        [nmax, y] = max(y,[],2) ;
    end
end