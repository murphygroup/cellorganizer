function [bestthresh] = ml_mlpthreshtest(netout,class,thresholds,onlyone)
%  ML_MLPTHRESHTEST - Threshold output of mlp for classification
%
%  [BESTTHRESH] = ML_MLPTHRESHTEST(NETOUT,CLASS,THRESHOLDS,ONLYONE)
%
%    Outputs:
%     BESTTHRESH - a vector of the best thresholds for all nets in netout
%                   where best = max(accuracy.^2 + recall.^2))
%
%    Inputs:
%     NETOUT - cell array of neural network outputs (rows=instances)
%     CLASS - cell array of one-of-N arrays indicating the true class of the 
%              outputs in NETOUT
%     THRESHOLD - vector of values to test as a threshold on the 
%              network outputs
%     ONLYONE - boolean value indicating whether the processing should 
%                consider a set of outputs to be unknown if more than one
%                of them is above the threshold.
%
%
%    M. Boland - 10 Apr 1999
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

% $Id: ml_mlpthreshtest.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

if(~iscell(netout) | ~iscell(class))
  error('NETOUT and CLASS must be cell arrays') ;
end

bestthresh = [] ;

%
% For each network...
for i=1:length(netout),

  accuracy   = [] ;
  recall     = [] ;

  %
  % Try each threshold with each network
  for j=1:length(thresholds),
    [acc rec] = ml_mlpthresh(netout{i},class{i},thresholds(j),onlyone) ;
    accuracy = [accuracy acc] ;
    recall   = [recall rec] ;
  end

  [maxval,maxidx] = max(accuracy.^2 + recall.^2) ;

  bestthresh = [bestthresh thresholds(maxidx)] ;

end
