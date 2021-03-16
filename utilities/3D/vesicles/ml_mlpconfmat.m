function [Cmat,Crate,missed]=ml_mlpconfmat(netout,class,threshold,onlyone)
% ML_MLPCONFMAT Compute a confusion matrix for output from a MLP.
%
% [CMAT,CRATE,MISSED]=ML_MLPCONFMAT(NETOUT,CLASS,THRESHOLD,ONLYONE)
%       Outputs:
%        CMAT - Confusion matrix in which the last column is 'unknown'
%        CRATE - five values: % of all samples classified correctly,
%                 % of attempts (total-unknown) classified correctly,
%                 the total number of samples,
%                 the number of samples classified correctly, and
%                 the number of samples classified as unknown.
%        MISSED - Structure with three components: a vector of indices
%                  indicating the misclassified samples, the true class
%                  of misclassified samples, and the assigned class
%                  of misclassified samples.
%
%       Inputs:
%        NETOUT - Matrix of network outputs
%        CLASS - Matrix of 1 of N classifications for the samples
%                 in NETOUT
%        THRESHOLD - Threshold to use when deciding which output 
%                     should be used for the classification
%        ONLYONE - boolean value indicating whether the processing should 
%                   consider a set of outputs to be unknown if more than one
%                   of them is above the threshold.
%
%
%	In the confusion matrix, the rows represent the true classes and the
%	columns the predicted classes.  
%
%  M. Boland - 13 Apr 1999
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

% $Id: ml_mlpconfmat.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

%
% numbers of outputs and classes
[numout numc]=size(netout);
[numout2 numc2]=size(class);

if numout~=numout2 | numc~=numc2
  error('Outputs and targets are different sizes')
end

if(numc==1)
  error('The classification array must be one-of-N encoding (i.e. each sample is specified by a row of N elements with a 1 in the position of the correct classification.') ;
end

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

% Compute 
numunknown = sum(threshclass==0) ;
correct = (trueclass==threshclass);
numcorrect=sum(sum(correct));
Crate=[numcorrect*100/numout numcorrect*100/(numout-numunknown) ...
       numout numcorrect numunknown];

%
% trueclass - true class number of each input
% imiss  - row index (into Y) of misclassified samples 
%
[val, imiss] = find(trueclass~=threshclass) ;
missed.index    = imiss ; 
missed.true     = trueclass(imiss) ; 
missed.assigned = threshclass(imiss) ;

%
% Generate the confusion matrix
Cmat=zeros(numc,numc+1);
for i=1:numc
  for j=1:numc
    Cmat(i,j) = sum((threshclass==j).*(trueclass==i));
    %
    % Extra class is unknown
    Cmat(i,numc+1) = sum((threshclass==0).*(trueclass==i)) ;
  end
end   
