function namelists = ml_classiflist(pred_y,y,setnames)
%ML_CLASSIFLIST Predicted classes in each true class.
%   NAMELISTS = ML_CLASSIFLIST(PRED_Y,Y,SETNAMES) returns a cell array of
%   the cell array of class names, which are predicted classes. PRED_Y is 
%   a vector of the predicted class labels and Y is a vector of true class
%   labels. SETNAMES is a cell array, in which each element is the name
%   of a class of the correponding label. For example, SETNAMES(i) is the
%   name of class i.

%   25-Jun-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University
 
if nargin < 3
    error('At least 3 arguments are required')
end

nclass=max(y);

for i=1:nclass
    pred_yi=pred_y(y==i);
    namelists{i}=setnames(pred_yi);
end
