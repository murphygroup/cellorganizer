function class = tz_lssvmclassify(sample,training,group,kernel)
%TZ_LSSVMCLASSIFY Least square SVM classification.
%   CLASS = TZ_LSSVMCLASSIFY(SAMPLE,TRAINING,GROUP,KERNEL) classify data
%   least square support vector machine (LSSVM). See CLASSIFY for
%   arguments SAMPLE, TRAINING, GROUP and CLASS. KERNEL is the kenel for
%   the surpport vector machine.

%   ??-???-???? Initial write T. Zhao
%   22-MAY-2004 ModifiedT. Zhao
%       - add kernel parameter
%   28-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('At least 4 arguments are required')
end

gam=10;
sig2=2;
type='classification';
model={training,group,type,gam,sig2,kernel,'preprocess'};

[alpha,b] = trainlssvm(model);
group=simlssvm(model,{alpha,b},sample);