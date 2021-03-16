function fullname = tz_getmodelfullname(modelname)
%ML_GETMODELFULLNAME Get full name of a regression model name
%   FULLNAME = ML_GETMODELFULLNAME(MODELNAME) returns a string representing
%   the name of the name abbreviation MODELNAME:
%   'bpnn' - neural network
%   'svm' - support vector machine
%   'lda - linear discriminant analysis
%   
%   See also

%   25-Jun-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required')
end

switch(modelname)
case 'bpnn'
    fullname='neural network';
case 'svm'
    fullname='support vector machine';
case 'lda'
    fullname='linear discriminant analysis';
otherwise
    fullname=modelname;
end