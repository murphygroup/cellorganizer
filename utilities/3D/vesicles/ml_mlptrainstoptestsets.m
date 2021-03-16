function [sets] = ml_mlptrainstoptestsets(nets,netin,trainidx,testidx,numsets,setsize,threshold,onlyone)
%  ML_MLPTRAINSTOPTESTSETS - Classify sets of instances by voting (plurality wins)
%
%  [SETS] = ML_MLPTRAINSTOPTESTSETS(NETS,NETIN,TRAINIDX,TESTIDX,NUMSETS,...
%            SETSIZE,THRESHOLD,ONLYONE) 
%
%    Outputs:
%     SETS - cell array of confusion matrices where the number of 
%             sets classified into each class is summed.  One entry
%             in the cell array per network.
%  
%    Inputs:
%     NETS - cell array of trained networks
%     NETIN - cell array of raw network inputs (NOT scaled mean=0, sdev=1)
%     TRAINIDX - cell array of indices used to generate training data
%     TESTIDX - cell array of indices used to generate test (stop) data
%     NUMSETS - number of sets to classify using each network
%     SETSIZE - size of each set
%     THRESHOLD - vector of values to test as a threshold on the 
%                  network outputs
%     ONLYONE - boolean value indicating whether the processing should 
%                consider a set of outputs to be unknown if more than one
%                of them is above the threshold.
%
%
%    M. Boland - 12 Apr 1999
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

% $Id: ml_mlptrainstoptestsets.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

if(~iscell(nets) | ~iscell(netin) | ~iscell(trainidx) | ~iscell(testidx))
  error('NETS, IDX, and NETIN must be cell arrays') ;
end

if (threshold ~= 0)
  if (length(threshold) ~= length(nets))
    error('The number of thresholds does not equal the number of networks.') ;
  end
else
  threshold = zeros(1,length(nets)) ;
end

%
% Reset the random number generator
rand('state',0) ;

%
% Number of classes
numc=nets{1}.nout ;

%
% Initialize the returned cell array
for i=1:length(nets)
  sets{i} = [] ;
end

for i=1:length(nets)
  trainall=[] ;
  %
  % Regenerate the training data for this network -- needed to normalize
  for j=1:numc
    trainall = [trainall ; netin{j}(trainidx{j}{i},:)] ;
  end

  for j=1:numc
    %
    % Normalize the test data for this class
    [trainnorm,testnorm] = ml_featurenorm(trainall,...
                                          netin{j}(testidx{j}{i},:)) ;
    sets{i} = [sets{i} ; ml_mlpsets(nets{i},testnorm,numsets,...
                                    setsize,threshold(i),onlyone)] ;
  end
end
