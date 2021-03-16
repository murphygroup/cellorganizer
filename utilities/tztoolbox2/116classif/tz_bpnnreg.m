function regmodel = tz_bpnnreg(x,y,t)
%TZ_BPNNREG Obsolete.

%function regmodel = tz_bpnnreg(x,y,t)
%OVERVIEW
%   regression by neural network
%PARAMETERS
%   x - variables
%   y - targets
%   t - regression parameters
%       t.prep: 0, normalize testing data based on all training data
%               1, normalize testing data based on training data without stop set
%RETURN
%   regmodel - trained model
%DESCRIPTION
%   
%HISTORY
%   16-May-2005 Initial write TINGZ
%SEE ALSO
%   tz_evalbpnnreg

error(tz_genmsg('of','tz_bpnnreg','ml_bpnnreg');

if ~exist('t','var')
    t.norm=1;
end

if ~isfield(t,'norm')
    t.norm=1;
end

if ~isfield(t,'hidden')
    t.hidden = 20;
end

if ~isfield(t,'epochs')
    t.epochs=300;
end

if ~isfield(t,'bp')
    t.bp=1000;
end

if ~isfield(t,'stop')
    t.stop=1;
end

if t.stop==1
    t.norm=2;
end

trainset=[];
stopset=[];
rlabel=[];
slabel=[];
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
    nsample=size(x,1);
    
    if ~isfield(t,'randtrainsel')
        t.randtrainsel=[];
    end
    if isempty(t.randtrainsel)
        t.randtrainsel=randperm(nsample);
    end
    
    stopsize=round(nsample/3);
    stopset=x(t.randtrainsel(1:stopsize),:);
    trainset=x(t.randtrainsel(stopsize+1:end),:);
    tmpytr=y(t.randtrainsel(stopsize+1:end),:);
    tmpyst=y(t.randtrainsel(1:stopsize),:);
    
    tmptrainset=trainset;
    
    switch(t.norm)
    case 1
        [trainset,zmean,zsdev]=tz_zscore(trainset);
        stopset=tz_zscore(stopset,zmean,zsdev);
        prep.zscore=1;
        prep.zmean=zmean;
        prep.zsdev=zsdev;
    case 2
        [trainset,zmean,zsdev]=tz_zscore(tmptrainset);
        stopset=tz_zscore(stopset,zmean,zsdev);
        prep.zscore=1;
        prep.zmean=zmean;
        prep.zsdev=zsdev;
    end
    
    % Train the network using the training and stop data
    [trainnetout stopnetout imin net] = ...
        tz_mlptraintest(trainset, tmpytr,stopset, tmpyst, t.hidden, t.epochs,t.bp) ;
else
    net = mlp(size(x,2), t.hidden, size(y,2), 'logistic') ;
    
    roptions = zeros(1,18) ;
    roptions(1) = 1 ;   % Output sse values
    %roptions(1) = -1 ;  % Output nothing 
    roptions(14) = t.epochs ;  % Number of epochs (train one epoch at a time)
    roptions(17) = 0.9 ;
    roptions(18) = 0.001 ;
       
    [net, roptions] = netopt(net, roptions, x, y, 'graddesc') ;     
end

regmodel=struct('modelname','bpnn','modeltype',...
    'network','trained',net,'t',t,'prep',prep);

%
% Summarize network performance
% y = mlpfwd(net, sample);
