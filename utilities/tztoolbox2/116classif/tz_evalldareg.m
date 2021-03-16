function y = tz_evalldareg(x,regmodel)
%TZ_EVALLDAREG Obsolete.

%function y = tz_evalldareg(x,regmodel)
%OVERVIEW
%   Evaluate lda regression
%PARAMETERS
%   x - evaluate lda regerssion
%   regmodel - trained lda regression model
%RETURN
%   y - predicted value
%DESCRIPTION
%   Only for categorized data. 
%HISTORY
%   29-May-2005 Initial write TINGZ
%SEE ALSO
%   tz_ldareg

error(ml_genmsg('of','tz_evalldareg','ml_evalldareg'));

if isfield(regmodel,'prep')
    if ~isempty(regmodel.prep.featidx)
        x=x(:,regmodel.prep.featidx);
    end
end

%number of groups
ngroup=length(regmodel.trained.rs);

%size of testing data
nx=size(x,1);

d = zeros(sr,ngroup);

for k = 1:ngroup
    meanx = regmodel.trained.means{k}(ones(nx,1),:);
    rinv = regmodel.trained.rs{k}'\(x-meanx)';
    d(:,k) = sum(rinv.*rinv)'*(regmodel.trained.rx(k)-1);
end
[mind, y] = min(d,[],2);
y=y';
