function [sorted, indices, distances, numpc, pcntexplained] = ...
    ml_TypICfm(AllFeatures, featselection, numpc, pcntexplained, robustcovar);
% [SCORES,INDICES,DISTANCES, NUMPC, PCNTEXPLAINED] =
%      ML_TYPICFM(ALLFEATURES,FEATSELECTION,NUMPC,PCNTEXPLAINED,ROBUSTCOVAR);
% ALLFEATURES   is the input feature matrix (columns=features, rows=images)
% FEATSELECTION identifies the method to be used for selecting a feature subset
%                'none'     - use all features calculated
%                'princomp' - calculate principal components and use in
%                              place of the features
% NUMPC         specifies how many principal components to use. (-1
%                to use PCNTEXPLAINED)
% PCNTEXPLAINED specifies the minimum percentage of variance to account for
% ROBUSTCOVAR   identified whether robust estimation of the
%                covariance matrix is to be used ('true' or 'false') outputs
% SCORES        scores of the images to rank    
% INDICES       contains the rankings of each image, with the most
%                typical image having the first (lowest) rank
% DISTANCES     contains the Mahalanobis distances of each image from the mean

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

% Modified Jennifer Lin from Kai's version of TypIXfm.m which was
%  orignially from Dr. Murphy's copy.
%
% July 2002 - EJSR - Formatting, changed variance sum loop
%                     fixed summing bug (numvars starts at 1)

[numfiles, numfeat] = size(AllFeatures);
numpc = floor(numpc);

switch featselection
 case 'none'
    featurematrix = double(AllFeatures);
    numvars = numfeat;
    pcntexplained = 100;
 case 'princomp'
    % Convert the features to zscores ((x-xmean)/std) so that the
    %  principal components are not dominated by features with large
    %  absolute values
    zfv = zscore(double(AllFeatures));
    % Now calculate principal components.  "pcs" will contain the pcs
    % themselves and variances will contain the amount of variance
    % explained by each pc
    [weights, pcs, variances, t2] = princomp(zfv);
    pcntvar = 100 * variances / sum(variances);
    if numpc<=0
	numvars = 0;
	currentTotal = 0;
	while(currentTotal < pcntexplained)
	    numvars = numvars +1;
	    currentTotal = currentTotal + pcntvar(numvars);
	end
	pcntexplained = currentTotal;
    else 
	numvars = numpc;
	if (numvars > size(AllFeatures, 2))
	    error(['Number of principal components requested is greater' ...
		   ' than the number of features calculated.']);
	end
	pcntexplained = sum(pcntvar(1:numvars));
    end
    % Make the feature matrix contain just the desired pcs
    featurematrix = pcs(:,1:numvars);
 otherwise
    error ('TypIC 2.0 - Invalid feature selection method specified!');
end


switch robustcovar
 case 'false'
    [distances, redcov, errortype] = ml_mahal(featurematrix, featurematrix);
 case 'true'
    [o1, o2, distances, o4, o5, o6] = ml_multout(featurematrix);
end

[sorted,indices] = sort(distances);
numpc = numvars;


 
