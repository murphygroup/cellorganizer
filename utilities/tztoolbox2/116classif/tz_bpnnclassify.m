function class = tz_bpnnclassify(sample,training,group,t)
%TZ_BPNNCLASSIFY Back propagation neural network (BPNN) classification.
%   CLASS = TZ_BPNNCLASSIFY(SAMPLE,TRAINING,GROUP) classifies the
%   data SAMPLE by BPNN. 1/3 of training data are used for stop set.
%   See CLASSIFY for more details about SAMPLE, TRAINING, GROUP and 
%   CLASS.
%   CLASS = TZ_BPNNCLASSIFY(SAMPLE,TRAINING,GROUP,T) lets users 
%   specify how to preprocess the data:
%       T = 0 or [] - z-score based on all training set
%       T = 1       - z-score without stop set
%       T = 2       - SDA and z-score without stop set

%   29-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('3 or 4 arguments are required')
end

if ~exist('t','var')
    t=0;
end

if isempty(t)
    t=0;
end

hidden = 20;
epochs=300 ;

caclass=ml_findclass(group);
nclass=length(caclass);
trainset=[];
stopset=[];
rlabel=[];
slabel=[];

if t==0
    [training,sample]=ml_featurenorm(training,sample);
end

for i=1:nclass
    ntrain=length(caclass{i});
    trainsize=round(ntrain*2/3);
    stopsize=round(ntrain-trainsize);
    trainset=[trainset;training(caclass{i}(1:trainsize,end),:)];
    rlabel=[rlabel;group(caclass{i}(1:trainsize,end),:)];
    stopset=[stopset;training(caclass{i}((trainsize+1):end,end),:)];
    slabel=[slabel;group(caclass{i}((trainsize+1):end,end),:)];
end

%remove constant features
% ntotal=size(trainset,2);
% fullidx=1:ntotal;
% constidx=tz_constidx(trainset);
% if ~isempty(constidx)
%     trainset(:,constidx)=[];
%     stopset(:,constidx)=[];
%     sample(:,constidx)=[];
%     fullidx(constidx)=[];
%     disp(['constant features removed:' num2str(constidx)]);
% end

tmptrainset=trainset;

switch(t)
case 1
    [trainset,stopset]=ml_featurenorm(tmptrainset,stopset);
    [trainset,sample]=ml_featurenorm(tmptrainset,sample);
case {2,3}
    ntotal=size(trainset,2);
    fullidx=1:ntotal;
    
    constidx=ml_constidx(trainset);
    
    if ~isempty(constidx)
        trainset(:,constidx)=[];
        stopset(:,constidx)=[];
        sample(:,constidx)=[];
        fullidx(constidx)=[];
        disp(['constant features removed:' num2str(constidx)]);
    end
    
    if t==2
        featsel=tz_sda(trainset,rlabel);
    else
        featsel=ml_sda(trainset,rlabel);
    end
    
    if isempty(featsel)
        warning('No good feature found. All features will be used');
        featsel=1:length(fullidx);
    end
    
    trainset=trainset(:,featsel);
    stopset=stopset(:,featsel);
    sample=sample(:,featsel);
    featsel=fullidx(featsel);
    
    save('/home/tingz/tmp/featsel.mat','featsel');
    tmptrainset=trainset;
    [trainset,stopset]=ml_featurenorm(tmptrainset,stopset);
    [trainset,sample]=ml_featurenorm(tmptrainset,sample);
end

% Train the network using the training and stop data
rpost=ml_label2post(rlabel);
spost=ml_label2post(slabel);

[trainnetout stopnetout imin net] = ...
    ml_mlptraintest(trainset, rpost.*0.8+0.1,stopset, spost.*0.8+0.1, hidden, epochs) ;

%
% Summarize network performance
netout = mlpfwd(net, sample);

% Find the largest output for each instance
[nmax, class] = max(netout,[],2) ;
