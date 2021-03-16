function s=tz_summaryobjclass(nfold,clstmeth,clustk,clstsda,classifmeth1, ...
    t1,classifmeth2,preprocess,t2,t3,featset,option,loadpath, ...
    savepath,other,merge)
%TZ_SUMMARYOBJCLASS Summary of object-level classification.
%   S = TZ_SUMMARYOBJCLASS(NFOLD,CLSTMETH,CLUSTK,CLSTSDA,CLASSIFMETH1, ...
%   T1,CLASSIFMETH2,PREPROCESS,T2,T3,FEATSET,OPTION,LOADPATH, ...
%   SAVEPATH,OTHER,MERGE) summarizes object-level classification
%   results. See TZ_TESTCV for details of the parameters.

%   ??-???-????  Initial write T. Zhao
%   05-NOV-2004  Modified T. Zhao
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 16
    error('Exactly 16 arguments are required')
end

clustername=[clstmeth num2str(clstsda) option other];
% if isempty(t3)
    setname=['cluster' clstmeth num2str(clstsda) '1st' classifmeth1 num2str(t1) '2nd' classifmeth2 num2str(preprocess) num2str(t2) featset option other];
% else
%     
%    setname=['cluster' clstmeth num2str(clstsda) '1st' classifmeth1 num2str(t1) '2nd' classifmeth2 num2str(preprocess) num2str(t2) 't2' num2str(t3) featset option other];
% end

kk=clustk;
for i=1:nfold
    clusterfilename=[num2str(kk) clustername num2str(i) 'fold.mat'];
    loadfile=[loadpath '/' num2str(kk) clustername num2str(i) 'fold.mat']
    load(loadfile);
    trainsels(i,:)=trainsel';
    testsels(i,:)=testsel';
    if exist('clstfeatsel','var')
        clstfeatsels{i}=clstfeatsel;
    else
        clstfeatsels={};
    end
    if any(trainsel+testsel)~=1
        warning('fold error')
    end
    
    ncluster(i)=max(tz_post2label(trainpost));
end

if any(sum(testsels,1)~=1)
    warning('fold error')
end

load([savepath '/' 'results' num2str(kk) setname num2str(i) 'fold.mat'])
[avgcm,avgacc,ncm,kappa,ua,pa]=tz_calcavgcm(ncvcm,merge);
if ~exist('sda','var')
    s=struct('nfold',nfold,'trainsels',trainsels,'testsels',testsels,...
        'ncluster',ncluster,'avgcm',avgcm,'avgacc',avgacc,'ncm',ncm,...
        'kappa',kappa,'ua',ua,'pa',pa,'aics',aic,'clstfeatsels',{clstfeatsels});
else
    
    s=struct('nfold',nfold,'trainsels',trainsels,'testsels',testsels,...
        'ncluster',ncluster,'avgcm',avgcm,'avgacc',avgacc,'ncm',ncm,...
        'kappa',kappa,'ua',ua,'pa',pa,'aics',aic,'clstfeatsels',{clstfeatsels},'sda',{sda});
end