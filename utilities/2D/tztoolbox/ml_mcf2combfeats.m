function [combfeats,combclass,combcellidx]=ml_mcf2combfeats(features)
%ML_MCF2COMBFEATS Convert MCF to combined feature matrix.
%   COMBFEATS = ML_MCF2COMBFEATS(FEATURES) returns the combined feature
%   matrix from the cell array of feature matrices FEATURES.
%   
%   [COMBFEATS,COMBCLASS,COMBCELLIDX] = ML_MCF2COMBFEATS(...) also returns
%   the combined class labels and cell indices.
%   
%   See also ML_COMBFEATS2MCF

%   AUG-07-2004 Initial write T. Zhao
%   31-OCT-2004 Modified T. Zhao
%       - add comments
%       - change function name tz_combinefeats_mcf --> tz_mcf2combfeats
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

nclass=length(features);
combfeats=[];
combclass=[];
combcellidx=[];

for i=1:nclass
    combfeats=[combfeats;features{i}];
    ncell=size(features{i},1);
    combclass=[combclass;zeros(ncell,1)+i];
    combcellidx=[combcellidx;(1:ncell)'];
end
