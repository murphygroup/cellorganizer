function y = tz_knn(xtrain, xtest, t, kmax)
%TZ_KNN	Creates a K-nearest-neighbour classifier.
%   Modified from KNN.
%	Y = KNN(XTRAIN, XTEST, T, KMAX) takes a matrix XTRAIN of input
%	vectors (one vector per row) with corresponding 1-of-N coded class
%	labels in the matrix T and uses the K-nearest-neighbour rule to
%	produce  a matrix Y of classifications for a new input matrix XTEST.
%	The nearest neighbours are determined using Euclidean distance. The
%	Kth column of the matrix Y contains the predicted class labels (as an
%	index 1..N, not as 1-of-N coding) using K nearest neighbours, where K
%	runs from 1 to KMAX.
%
%	See also
%	KMEANS
%

%	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)
%   Modified by tingz on Apr. 30, 2004
%    - Only do classification for kmax 

if nargin < 4
    error('At least 4 arguments are required')
end

ntrain = size(xtrain, 1);
% Check that dimensions of datasets are consistent
if size(xtrain, 2) ~= size(xtest, 2)
  error('Inconsistent number of variables in training and test data.')
end

if ntrain ~= size(t, 1)
  error('Inconsistent number of examples in input and target data.')
end

ntest = size(xtest, 1);		% Number of test vectors.
nclass = size(t, 2);		% Number of classes.

% Compute matrix of squared distances between input vectors from the training 
% and test sets.  The matrix distsq has dimensions (ntrain, ntest).

distsq = dist2(xtrain, xtest);

% Now sort the distances. This generates a matrix kind of the same 
% dimensions as distsq, in which each column gives the indices of the
% elements in the corresponding column of distsq in ascending order.

[vals, kind] = sort(distsq);
y = zeros(ntest, kmax);
sums = zeros(ntest, nclass);

for k=1:kmax    %original: k=1:kmax
  % We now look at the predictions made by the Kth nearest neighbours alone,
  % and represent this as a 1-of-N coded matrix, and then accumulate the 
  % predictions so far.

  sums = sums + t(kind(k,:),:);

  % Now pick out the majority vote.
  % To break ties in a random fashion, we add a small random component.
  % The vector INDEX of length NTEST contains the index of the winning
  % class. 

  [temp, index] = max(sums' + 0.5*rand(size(sums')));

  % Now load this into the Kth column of Y.
  y = index';  %original: y(:,k) = index';
end
