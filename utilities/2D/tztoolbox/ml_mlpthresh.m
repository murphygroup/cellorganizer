function [accuracy, recall] = ml_mlpthresh(netout,class,threshold,onlyone)
%  ML_MLPTHRESH - Threshold output of mlp for classification
%
%  [ACCURACY RECALL] = ML_MLPTHRESH(NETOUT,CLASS,THRESHOLD,ONLYONE)
%
%    Outputs:
%     ACCURACY - number of correct classifications / total classifications
%                 attempted (not including unknowns)
%     RECALL - number of correct classifications / total instances
%
%    Inputs:
%     NETOUT - output of neural network (rows=instances)
%     CLASS - one-of-N array indicating the true class of the 
%              outputs in NETOUT
%     THRESHOLD - value to use as a threshold on the network outputs
%     ONLYONE - boolean value indicating whether the processing should 
%                consider a set of outputs to be unknown if more than one
%                of them is above the threshold.
%
%
%    M. Boland - 01 Mar 1999
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

% $Id: ml_mlpthresh.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

%
% How many classes?
numc = size(netout,2) ;

%
% Convert the 1 of N format of class to single digits per instance
trueclass = [1:numc]*class' ;

%
% Find those outputs that are above the threshold
netthresh = netout .* (netout>=threshold) ;

%
% Find the largest output for each instance
[nmax, threshclass] = max(netthresh') ;

%
% Don't consider classifications where max=0
threshclass = threshclass .* (nmax>0) ;

%
% Identify those sets of outputs for which only one value is above threshold
if(onlyone)
  justone = (sum((netout>=threshold)') == 1) ;
  threshclass=justone.*threshclass ;
end

%
% Number of correctly classified samples
numcorrect = sum(threshclass == trueclass) ;

numinstances = size(netout,1) ;
numattempted = sum(threshclass>0) ;

if(numattempted>0),
  accuracy = numcorrect/numattempted ;
else
  accuracy=0 ;
end

if(numinstances>0),
  recall = numcorrect/numinstances ;
else
  recall = 0 ;
end
