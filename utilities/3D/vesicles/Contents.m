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

% Murphy Lab Classification
%
% General NN
%   ml_make_nn_classifier  - Trains and Runs Neural Network
%                            Classifier
%                            Currently returns confmats from ml_confmat
%
% NN implementation
%   ml_mlptrainstoptest    - Wrapper to Train network with n-fold
%                            cross-validation, using a stop set to
%                            complete training
%   ml_mlptraintest        - Trains network and tests it. The testing
%                            actually determines at what stage a
%                            minimum value is reached and retrains up
%                            to this stage (i.e test with stop set)
%   ml_mlptrain            - Trains network for each epoch
%   ml_featurenorm         - Normalizes features for training
%
% NN output
%   ml_confmat             - Analysis of Training and testing of a
%                            NN. Winner is the maximum output value.
%                            Returns the confusion matrix of a single
%                            network result. (Confusion matrix is [10x10])
%   ml_mlpconfmat          - Analysis of training and testing of a
%                            NN. Winner is the output which is greater
%                            than a given threshold. If 0
%                            outputs are above threshold,
%                            classification is judged as unknown
%                            (confusion matrix is [10x11])
%                            If more than one output is above
%                            threshold, either the max can be chosen,
%                            or the instance can be classified as unknown
%   ml_mlpconfmatall       - Analysis of training and testing of
%                            multiple NN (i.e permutations). Calls
%                            ml_mlpconfmat for each network.
%   ml_mlptrainstoptestsets- Analysis of training by testing groups
%                            (sets) of instances, classifying by
%                            plurality
%   ml_mlpclasssummary     - Taking all confusion matrices from
%                            ml_mlpconfmatall, a summary will be
%                            generated of the average confusion
%                            matrix, mean and variance of true
%                            classifer performance (for confusion
%                            matrix of all classes, and confusion
%                            matrix of all classes + unknown)
%   ml_mlpthresh           - Determine accuracy and recall at a given
%                            threshold. Choose a threshold for which
%                            accuracy^2+recall^2 is maximized.
%   ml_mlpsetsummary       - Taking all confusion matrices from 
%                            (as a cell arrary), a summary will be 
%                            generated of the average confusion
%                            matrix without unknowns, the average
%                            with unknowns and the classifier
%                            performance. Note, this is the same as
%                            ml_mlpclasssummary with different
%                            format of input and output
%
% NN classification
%   ml_mlpsets             - Run a randomized set of instances from a
%                            single class to determine what that class
%                            is. Classification is based on plurality
%                            (i.e. majority predicted class is the
%                            winner). The instance features MUST be
%                            normalized (mean of 0, variance of 1). It
%                            is also possible to test a single
%                            instance.
%
