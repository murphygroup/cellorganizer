function class = tz_svmclassify(sample,training,group,param)
%TZ_SVMCLASSIFY Support vector machine (SVM) classifier.
%   CLASS = TZ_SVMCLASSIFY(SAMPLE,TRAINING,GROUP,PARAM) classifies the
%   data SAMPLE by the support vector machine toolbox downloaded from
%   http://www.ece.osu.edu/~maj/osu_svm/ . See CLASSIFY for more details
%   about SAMPLE, TRAINING, GROUP and CLASS. PARAM contains SVM
%   parameters. See MEXSVMTRAIN for more details.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('Exactly 4 arguments are required')
end

[AlphaY, SVs, Bias, Parameters, nSV, nLabel] = ...
    mexSVMTrain(training', group',param);
[ClassRate, DecisionValue, Ns, ConfMatrix, PreLabels]= ...
    mexSVMClass(sample', ones(1,size(sample,1)), AlphaY, ...
    SVs, Bias, Parameters, nSV, nLabel);
class=PreLabels';