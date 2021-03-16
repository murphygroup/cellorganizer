function [rnetout,tnetout,imin,net]=ml_mlptraintest(ri,ro,ti,to,hidden,epochs,dp)
%  ML_MLPTRAINTEST - Train and test a multilayer perceptron
%  [RNETOUT,TNETOUT,IMIN,NET] = ML_MLPTRAIN(RI,RO,TI,TO,HIDDEN,DP)
%
%  rnetout - output from trained net for training samples
%  tnetout - output from trained net for test samples
%  imin    - epoch at which the sum of squared error was minimized for the 
%            test set.
%  net     - network after training
%
%  ri      - tRaining set inputs. samples are rows, features are columns
%  ro      - tRaining set outputs . samples are rows, net outputs are columns
%  ti, ti  - same for the Test data
%  hidden  - number of hidden nodes
%  epochs  - number of epochs to complete before checking for
%             a minimum in the SSE on the test (stop) data
%  dp      - 1/dp is the learning rate. default 1000
%  Default parameters:
%   momentum = 0.9

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

% $Id: ml_mlptraintest.m,v 1.4 2006/06/27 13:33:47 tingz Exp $

if ~exist('dp','var')
    dp=1000;
end

net = ml_mlp(size(ri,2), hidden, size(ro,2), 'logistic') ;

roptions = zeros(1,18) ;
roptions(1) = 1 ;   % Output sse values
%roptions(1) = -1 ;  % Output nothing 
roptions(14) = 1 ;  % Number of epochs (train one epoch at a time)
roptions(17) = 0.9 ;
roptions(18) = 0.001 ;

[net,rsse,tsse] = ml_mlptrain(net, roptions, ri, ro, ti, to, epochs) ;

%
% Round tsse to the nearest 1/dp to avoid long training sessions
%  with little progress
%
[minv,imin] = min(round(tsse*dp)/dp) ;
pass=1 ;

%
% Continue training until the minimum SSE occurs less than 90% 
%  of the way through the last pass.
%
while ((imin ./ pass) > 0.9*epochs)
  rssesave = rsse ;
  tssesave = tsse ;
  [net,rsse,tsse] = ml_mlptrain(net,roptions,ri,ro,ti,to,epochs) ;
  rsse = [rssesave rsse] ;
  tsse = [tssesave tsse] ;
  [minv,imin] = min(round(tsse*dp)/dp) ;
  %plot([1:length(rsse)], rsse, 'r-', [1:length(tsse)], tsse, 'b-');
  pass=pass+1 ;
  pass
  minv
  imin
end

%plot([1:length(rsse)], rsse, 'r-', [1:length(tsse)], tsse, 'b-') ;

net = ml_mlp(size(ri,2), hidden, size(ro,2), 'logistic') ;

[net,rsse,tsse] = ml_mlptrain(net, roptions, ri, ro, ti, to, imin) ;

rnetout = mlpfwd(net, ri) ;
tnetout = mlpfwd(net, ti) ;

