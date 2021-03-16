function regmodel = tz_svmreg(x,y,t)
%TZ_SVMREG Obsolete

%function regmodel = tz_svmreg(x,y,t)
%OVERVIEW
%   regression by svm
%PARAMETERS
%   x - variables nxm
%   y - targets nx1 or nxk (k classes)
%       e.g. y=[1 1 2 2 3 3]';
%           is the same as
%           y=[1 -1 -1;1 -1 -1;-1 1 -1;-1 1 -1;-1 -1 1;-1 -1 1];
%   t - regression parameters
%       t.norm: 0, normalize testing data based on all training data
%               1, normalize testing data based on training data without stop set
%RETURN
%   regmodel - trained model
%DESCRIPTION
%   Onlyt take categorized data
%HISTORY
%   15-June-2005 written by SHANNCC
%SEE ALSO
%   tz_evalsvmreg

error(tz_genmsg('of','tz_svmreg','ml_svmreg'));

if ~exist('t','var')
    t.norm=1;
end

if ~isfield(t,'norm')
    t.norm=1;
end

if ~isfield(t,'stop')
    t.stop=0;
end

prep.zmean=[];
prep.zsdev=[];
prep.zscore=0;
constidx=tz_constidx(x);
if isempty(constidx)
    prep.featidx=[];
else
    allidx=1:size(x,2);
    allidx(constidx)=[];
    prep.featidx=allidx;
    x(:,constidx)=[];
end

%normalize upon all training data
if t.norm==1
    [x,zmean,zsdev]=tz_zscore(x);
    prep.zscore=1;
    prep.zmean=zmean;
    prep.zsdev=zsdev;
end

if t.stop==1

else
    net = t.model_types;
    net = train(t.model_types, t.tutor, x, y, t.C_values, rbf(t.rbf_levels));    
end

regmodel=struct('modelname','svm','modeltype',...
    'svm','trained',net,'t',t,'prep',prep);

if size(y,2)==1
    regmodel.postp.ctg=1;
end