function [Cmat,Crate,missed]=ml_confmat(Y,T)
%ML_CONFMAT Compute a confusion matrix
%
% [CMAT, CRATE, MISSED] = ML_CONFMAT(Y, T)
%      Outputs:
%        CMAT   - Confusion matrix [NumClasses, NumClasses]
%        CRATE  - Classification performance [% correct, Num correct]
%        MISSED - Structure with three components: a vector of indices
%                  indicating the misclassified samples, the true class
%                  of misclassified samples, and the assigned class
%                  of misclassified samples
%      Inputs:
%        Y      - Single Network output [NumInstances x NumOutputs]
%        T      - Actual instance classes [NumInstances x NumOutputs]

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

%Based on:
%CONFMAT Compute a confusion matrix.
%
%	Description
%	[CMAT, CRATE,MISSED]=MB_CONFMAT(Y,T) computes the confusion matrix Cmat
%	and classification performance RATE for the predictions mat{y} compared
%	with the targets T.  The data is assumed to be in a 1-of-N encoding,
%	unless there is just one column, when it is assumed to be a 2 class
%	problem with a 0-1 encoding.  Each row of Y and T corresponds to a
%	single example.
%
%	In the confusion matrix, the rows represent the true classes and the
%	columns the predicted classes.  The vector RATE has two entries: the
%	percentage of correct classifications and the total number of correct
%	classifications.
%
%	See also
%	CONFFIG, DEMTRAIN
%
%	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)
%
% Modified 26 jan 1999 by M. Boland
%

% $Id: ml_confmat.m,v 1.3 2006/06/27 13:33:47 tingz Exp $


[n c]=size(Y);
[n2 c2]=size(T);

if n~=n2 | c~=c2
    error('Outputs and targets are different sizes')
end

if c > 1
    % Find the winning class assuming 1-of-N encoding
    [maximum Yclass] = max(Y', [], 1);
    
    TL=[1:c]*T';
else
    % Assume two classes with 0-1 encoding
    c = 2;
    class2 = find(T > 0.5);
    TL = ones(n, 1);
    TL(class2) = 2;
    class2 = find(Y > 0.5);
    Yclass = ones(n, 1);
    Yclass(class2) = 2;
end

% Compute 
correct = (Yclass==TL);
total=sum(sum(correct));
Crate=[total*100/n total];

%
% M. Boland -- not tested with 2 class input
%
% trueclass - true class number of each input
% imiss  - row index (into Y) of misclassified samples 
%
[maxc, trueclass] = max(T',[],1) ;
[val, imiss] = find(Yclass~=TL) ;
missed.index    = imiss ; 
missed.true     = trueclass(imiss) ; 
missed.assigned = Yclass(imiss) ;

Cmat=zeros(c,c);
for i=1:c
    for j=1:c
        Cmat(i,j) = sum((Yclass==j).*(TL==i));
    end
end   
