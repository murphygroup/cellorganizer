function [ncm,pcm] = ml_summarytestclassif(y,pred_y,tr_y)
%ML_SUMMARYTESTCLASSIF Summarize classification for testing data.
%   NCM = ML_SUMMARYTESTCLASSIF(Y,PRED_Y,TR_Y) returns a confusion matrix
%   with numbers by checking true class labels Y and predicted class
%   labels PRED_Y, which are labels in TR_Y. Y, PRED_Y and TR_Y must be
%   vectors.
%
%   [NCM,PCM] = ML_SUMMARYTESTCLASSIF(Y,PRED_Y,TR_Y) also returns
%   confusion matrix with percentages.
%
%   See also ML_SUMMARYCLASSIF

%   15-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University
  
if nargin < 3
    error('Exactly 3 arguments are required')
end

if all( size(y) > 1 )
    error('The 1st argument is not a vector.');
end

if all( size(pred_y) > 1 )
    error('The 2nd argument is not a vector.');
end

if all( size(tr_y) > 1 )
    error('The 3rd argument is not a vector.');
end

if length(y) ~= length(pred_y)
    error('The two arguments do not have the same length.');
end

ntestclass=max(y);
ntrainclass=max(tr_y);

for i=1:ntestclass
    for j=1:ntrainclass
        index=find(y==i);
        if isempty(index)
            ncm(i,j)=0;
        else
            ncm(i,j)=sum(double(pred_y(index)==j));
        end
    end
end

pcm=ml_normrow(ncm)*100;

