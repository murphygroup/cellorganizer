function [y,prob] = ml_evalreg(x,regmodel,postp)

%ML_EVALREG evaluates the trained regression model
%   Y=ML_EVALREG(X,REGMODEL) return the predicted target
%   value of X according REGMODEL, which is returned from
%   ML_REGRESS.
%   Y=ML_EVALREG(X,REGMODEL,POSTP) takes POSTP as a post
%   processing method for predicted data. If POSTP.ctg is 1, 
%   it will return labels. And this is the only method supported
%   by current version.
%
%   See also ml_regress, ml_trainclassif
  
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

%HISTORY
%   12-Jul-2005 Initial write TINGZ

if isfield(regmodel,'succ')
     if regmodel.succ==0
         warning('An invalid model is taken.');
         y=0;
         return;
     end
end

if nargin>2
    regmodel.postp=postp;
end

% if isfield(regmodel.t,'idx_sda')
%     x = x(:,regmodel.t.idx_sda);
% end

switch regmodel.modelname
    case 'bpnn'
        if nargout==2
            [y, prob] = ml_evalbpnnreg(x,regmodel);
        else
            y = ml_evalbpnnreg(x,regmodel);
        end
    case 'svm'  % prob is raw score, not the real prob
        if nargout==2
            [y, prob] = ml_evalsvmreg(x,regmodel);
        else
            y = ml_evalsvmreg(x,regmodel);
        end
    case 'lda' % prob is d, still wait for validation
        [y, d] = ml_evalldareg(x,regmodel);
        if nargout==2
            unnormprob = exp(-d);
            prob = unnormprob ./ repmat(sum(unnormprob, 2), [1 size(unnormprob,2)]);
        end
    case 'libsvm'
        if nargout==2
            [y, prob] = ml_evallibsvmreg(x,regmodel);
        else
            y = ml_evallibsvmreg(x,regmodel);
        end
end
