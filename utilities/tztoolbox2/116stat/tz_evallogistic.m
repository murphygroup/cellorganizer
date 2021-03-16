function p=tz_evallogistic(x,para)
%TZ_EVALLOGISTIC Obsolete. See ML_EVALOGISTIC.
%   P = TZ_EVALLOGISTIC(X,PARA)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function p=tz_evallogistic(x,para)
%
%OVERVIEW:
%   get probabilities of the logistic model
%PARAMETERS:
%   x - data, each column is a variable
%   para - parameters
%RETURN:
%   p - probability
%DESCRIPTION:
%
%HISTORY:
%   ??-OCT-2004 Initial write TINGZ
%   31-OCT-2004 Modified TINGZ
%       - add comments

error(tz_genmsg('of','tz_evallogistic','ml_evallogistic'));

y=x*para;

p=1./(1+exp(-y));
