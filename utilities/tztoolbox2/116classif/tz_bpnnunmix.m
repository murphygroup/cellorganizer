function membts = tz_bpnnunmix(sample,training,membtr,t)

%ML_BPNNUNMIX Unmixing by BPNN.
%   MEMBTS = ML_BPNNUNMIX(SAMPLE,TRAINING,MEMBTR,T) returns the proportions
%   of testing samples SAMPLE in all components. The number of components
%   is the same number of MEMBTS and MEMBTR.
%   TRAINING is the training data, which must have the same
%   number of rows as MEMBTR. The values in MEMBTR must be between [0,1]
%   the sum of each row must be 1. T is the training parameters, see
%   ML_REGRESS for more details.
%   
%   See also ML_REGRESS

%   13-May-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University  

if size(membtr,2)<=1
    error('the number of columns of membtr must be greater than 1');
end

if any(membtr<0) | any(membtr>1) | any(abs(sum(membtr,2)-1)>1e-5)
    error('invalid proportions');
end

ytn=membtr(:,1:end-1);
regmodel = tz_bpnnreg(training,ytn,t);
ytt=tz_evalbpnnreg(sample,regmodel);

membts=[ytt,1-sum(ytt,2)];
membts(membts<0)=0;
membts=tz_normobjcom(membts,0)