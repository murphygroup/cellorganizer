function [comb,combclass,combcellidx] = ml_mc2comb(mc)
%ML_MC2COMB
%   COMB = ML_MC2COMB(MC)
%   
%   [COMB,COMBCLASS,COMBCELLIDX] = ML_MC2COMB(...)
%   
%   See also

%   27-Oct-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 1
    error('Exactly 1 argument is required');
end

nclass=length(mc);
comb={};
combclass=[];
combcellidx=[];

for i=1:nclass
    comb={comb{:} mc{i}{:}};
    ncell=length(mc{i});
    combclass=[combclass;zeros(ncell,1)+i];
    combcellidx=[combcellidx;(1:ncell)'];
end
