function [net, rsse, tsse] = ml_mlptrain(net, roptions, ri, ro, ti, to, epochs)
%  ML_MLPTRAIN(NET,ROPTIONS,RI,RO,TI,TO,EPOCHS)
%  [NET,RSSE,TSSE] = ML_MLPTRAIN(NET,ROPTIONS,RI,RO,TI,TO,EPOCHS)
%
%    Outputs:
%     NET - neural network after training
%     RSSE - vector of sum of squared error for the training data
%     TSSE - vector of sum of squared error for the test data
%
%    Inputs:
%     NET - neural network before training
%     ROPTIONS - 18 element vector of net options
%     RI - training data inputs
%     RO - training data outputs
%     TI - test data inputs
%     TO - test data outputs
%     EPOCHS - number of training epochs
%
%    To generate NET:
%     net=mlp(nin, nhidden, nout, func), where func is 'linear',
%         'logistic', or 'softmax'.
%
%    To generate ROPTIONS:
%     roptions = zeros(1,18) ; 
%     roptions(14) = 1 ;     % Number of iterations
%     roptions(1) = 1 ;      % Does not work unless this option is set hmmm...
%     roptions(17) = 0.9 ;   % Momentum
%     roptions(18) = 0.001 ; % Learning Rate
%
%    After training:
%     rnetout = mlpfwd(net, ri) ;
%     tnetout = mlpfwd(net, ti) ;
%     [rcmat, rcrate, rmissed] = ml_confmat(rnetout, rclass) ;
%     [tcmat, tcrate, tmissed] = ml_confmat(tnetout, tclass) ;
%
%    M. Boland - 16 Feb 1999
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

% $Id: ml_mlptrain.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

tsse = [] ;
rsse = [] ;

for i = 1:epochs,
	[net, roptions] = netopt(net, roptions, ri, ro, 'graddesc') ;
	rsse = [rsse roptions(1,8)] ;
	tsse = [tsse mlperr(net, ti, to)] ;
end

