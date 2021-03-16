function [nn, trainnetout, stopnetout, testnetout, trainclasses, ...
	 stopclasses, testclasses, trainidx, stopidx, testidx, imin] = ...
    ml_make_nn_classifier( features, subset_index, TrainSetSize, ...
				StopSetSize, HiddenLayerSize, ...
				EpochsAtOnce, NoOfNets)
% ML_MAKE_NN_CLASSIFIER    Train and test a Netlab backprop net
%    [NN, trainnetout, stopnetout, testnetout, trainclasses,...
%     stopclasses, testclasses, trainidx, stopidx, testidx, imin] = ...
%                   ML_MAKE_CLASSIFIER( FEATURES,
%                   SUBSET_INDEX, TRAINSETSIZE, STOPSETSIZE,
%                   HIDDENLAYERSIZE, EPOCHSATONCE, NOOFNETS)
%

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




if (nargin < 7)
    NoOfNets = 10; %a.k.a. Permutations
end
if (nargin < 6)
    EpochsAtOnce = 1;
end
if (nargin < 5)
    HiddenLayerSize = 20;
end
if (nargin < 4)
   StopSetSize = 20;
end
if (nargin < 3)
    TrainSetSize = 40;
end
if (nargin < 2)
    subset_index = [1:size(features{1},2)];
end
if (nargin < 1)
    disp('Must input features');
    return;
end



% make arrays for Desired Outputs and Class Assignments
for c = 1 : length( features)
  features{c} = double(features{c}(:,subset_index));
  %class_assignments{c} = zeros(size(features{c},1),10);
  class_assignments{c} = zeros(size(features{c},1), length(features));
  class_assignments{c}(:,c) = 1;
  desired_outputs{c} = class_assignments{c} * 0.8 + 0.1;
end


[trainnetout,stopnetout,testnetout,trainclasses,stopclasses,testclasses, ...
trainidx,stopidx,testidx,imin,nn] = ml_mlptrainstoptest( ...
    features,desired_outputs,class_assignments, ...
    HiddenLayerSize, ...
    TrainSetSize,StopSetSize, ...
    EpochsAtOnce,NoOfNets);

% Make confusion matrix for test set
crates=[];
for i=1:length(testnetout)
 [CMAT, CRATE,MISSED]=ml_confmat(testnetout{i}, testclasses{i});
 cmat(:,:,i) = CMAT;
 crates = [crates CRATE(1)];
end
avg_cmat = mean( cmat, 3);
perc_mat = 100*avg_cmat./repmat(sum(avg_cmat,2),[1 size(avg_cmat,2)]);
confmat = perc_mat;
confmats = cmat;
