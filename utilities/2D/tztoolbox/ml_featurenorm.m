function [train_norm, test_norm] = ml_featurenorm(train, test)
% ML_FEATURENORM - Normalize training and test data
%
% [TRAIN_NORM, TEST_NORM] = ML_FEATURENORM(TRAIN, TEST)
%    Normalizes the data in train to have mean=0 and variance=1.
%    The values used to normalize train are also used to normalize test.
%
% M. Boland - 18 Jan 1999

% Copyright (C) 2006  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

% $Id: ml_featurenorm.m,v 1.3 2006/06/27 13:33:47 tingz Exp $

%
% R - number of training samples
% T - number of test samples
%
R = size(train,1) ;
T = size(test,1) ;

%
% Avoid division by 0 by mapping var~0 to var=1
%  NOTE: empirically not necessary -- added after processing
%         an incorrect data set (i.e. a single class)
train_sdev = sqrt(var(train)) ;

%check constant index
constidx=ml_constidx(train);
if ~isempty(constidx)
    train_sdev(constidx)=1;
end

% train_sdev = train_sdev + (train_sdev < 1e-10 ) ;
%
% Why was this here?
%test_sdev = sqrt(var(test)) ;
% test_sdev = test_sdev + (test_sdev < 1e-10) ;

train_norm = (train-(ones(R,1)*mean(train))) ./ (ones(R,1)*train_sdev);
if (T > 0)
  test_norm = (test-(ones(T,1)*mean(train))) ./ (ones(T,1)*train_sdev);
else
  test_norm = [] ;
end

