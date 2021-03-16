function features=ml_combfeats2mcf(combfeats,combclass)
%ML_COMBFEATS2MCF Change feature structure from combined to MCF.
%   FEATURES = ML_COMBFEATS2MCF(COMBFEATS,COMBCLASS) returns MCF from the
%   combined feature matrix COMBFEATS and class labels COMBCLASS, which
%   is a column vector. COMBFEATS and COMBCLASS must have the same number
%   of rows.
%   
%   See also ML_MCF2COMBFEATS

%   ??-???-???? Initial write T. Zhao
%   31-OCT-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% Copyright (C) 2007  Murphy Lab
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

if nargin < 2
    error('Exactly 2 arguments are required')
end

features=ml_findclass([combfeats,combclass]);
nclass=length(features);

for i=1:nclass
    features{i}=features{i}(:,1:end-2);
end