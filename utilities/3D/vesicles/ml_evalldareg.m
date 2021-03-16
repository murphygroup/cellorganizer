function [y,d] = ml_evalldareg(x,regmodel)

%ML_EVALLDAREG classifies the input data by a trained model
%   Y=ML_EVALLDAREG(X,REGMODEL) returns the group labels of X. For example,
%   if there are n classes, Y will be a vector of values from 1 to n
%   REGMODEL is a structure returned from ML_LDAREG.
%
%   [Y,D]=ML_EVALLDAREG(...) also returns the log of minus postirior l
%   ikelihood.
%
%   See also ML_LDAREG

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

%   29-May-2005 Initial write TINGZ

if isfield(regmodel,'prep')
    if ~isempty(regmodel.prep.featidx)
        x=x(:,regmodel.prep.featidx);
    end
end

%number of groups
ngroup=length(regmodel.trained.rs);

%size of testing data
nx=size(x,1);

d = zeros(nx,ngroup);

if isfield(regmodel.trained,'prior')
    for k = 1:ngroup
        meanx = regmodel.trained.means{k}(ones(nx,1),:);
        rinv = (x-meanx)/regmodel.trained.rs{k};
        d(:,k) = log(regmodel.trained.prior(k)) - ...
            .5*(sum(rinv.*rinv,2)+regmodel.trained.logsigma(k));
    end
    [maxd, y] = max(d,[],2);
    d = -d;  % negative log likelihood
else %old version
    for k = 1:ngroup
        meanx = regmodel.trained.means{k}(ones(nx,1),:);
        rinv = regmodel.trained.rs{k}'\(x-meanx)';
        d(:,k) = sum(rinv.*rinv)'*(regmodel.trained.rx(k)-1);
    end
    [mind, y] = min(d,[],2);
end


