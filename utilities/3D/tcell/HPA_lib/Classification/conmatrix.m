function [C,U,D,W] = conmatrix( truelabels,predlabels, Nclasses)
% [CMAT ACC WCMAT WACC] = CONFMAT( TRUELABELS, PREDLABELS)
%    calculates the confusion matrix. CMAT is the matrix, 
%    ACC is the classification accuracy, WCMAT is the 
%    weighted confusion matrix, and WACC is the weighted 
%    classification accuracy.
% 
% Last modified Justin Newberg 23 October 2007

if ~exist('Nclasses','var')
    Nclasses = 0;
end

MI = 1;
MA = max(max(truelabels),Nclasses);
C = zeros(MA,MA);
for i=1:MA
	idx = find( truelabels==i);
	[c b] = hist( predlabels( idx), [MI:1:MA]);
	C(i,:) = c;
end

U = sum( reshape( (eye(MA,MA).*C), [MA*MA 1])) / sum(C(:));

D = C./repmat( sum(C,2), [1 MA]);

D(isnan(D))=0;
W = sum(sum(eye(MA).*D)) / (MA-sum(isnan(sum(D,2))));
